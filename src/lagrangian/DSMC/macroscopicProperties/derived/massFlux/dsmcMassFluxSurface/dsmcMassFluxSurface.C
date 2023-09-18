/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dsmcMassFluxSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dsmcMassFluxSurface, 0);
addToRunTimeSelectionTable(dsmcField, dsmcMassFluxSurface, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcMassFluxSurface::dsmcMassFluxSurface
(
    const Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.get<word>("field")),
    regionId_(-1),
    sampleInterval_(propsDict_.getOrDefault<label>("sampleInterval", 1)),
    sampleCounter_(0),
    faceZoneName_(propsDict_.get<word>("faceZone")),
    zoneSurfaceArea_(0.0),
    typeIds_(),
    fluxDirection_(propsDict_.get<vector>("fluxDirection")),
    molsZone_(0.0),
    massZone_(0.0),
    averagingCounter_(0.0),
    timeIndex_(0),
    molFluxZone_(1, 0.0),
    massFluxZone_(1, 0.0),
    massFlowZone_(1, 0.0),
    averagingAcrossManyRuns_
    (
        propsDict_.getOrDefault<bool>("averagingAcrossManyRuns", false)
    )
{
    // choose molecule ids to sample
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Read stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        Info << nl << "Averaging across many runs initiated." << nl << endl;

        read();
    }


    // select face zone

    const faceZoneMesh& faceZones = mesh_.faceZones();

    regionId_ = faceZones.findZoneID(faceZoneName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << faceZoneName_ << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    fluxDirection_.normalise();

    // find total surface area
    const labelList& faces = faceZones[regionId_];

    if (Pstream::parRun())
    {
        DynamicList<label> processorFaces(0);

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                for (label p = 0; p < patch.size(); ++p)
                {
                    label patchFaceI = p + patch.start();
                    label faceId = faces.find(patchFaceI);

                    if (faceId != -1)
                    {
                        processorFaces.append(patchFaceI);
                    }
                }
            }
        }

        processorFaces.shrink();

        label nInternalFaces = faces.size() - processorFaces.size();

        List<label> internalFaces(nInternalFaces, 0);

        label counter = 0;

        for (const label faceI : faces)
        {
            if (processorFaces.find(faceI) == -1)
            {
                internalFaces[counter] = faceI;
                ++counter;
            }
        }

        for (const label faceI : internalFaces)
        {
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }

        for (const label faceI : processorFaces)
        {
            zoneSurfaceArea_ += 0.5*mag(mesh_.faceAreas()[faceI]);
        }

        if (Pstream::parRun())
        {
            for (int p = 0; p < Pstream::nProcs(); ++p)
            {
                if (p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << zoneSurfaceArea_;
                    }
                }
            }

            // receiving
            for (int p = 0; p < Pstream::nProcs(); ++p)
            {
                if (p != Pstream::myProcNo())
                {
                    scalar zoneSurfaceAreaProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> zoneSurfaceAreaProc;
                    }

                    zoneSurfaceArea_ += zoneSurfaceAreaProc;
                }
            }
        }
    }
    else
    {
        forAll(faces, f)
        {
            const label faceI = faces[f];

            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }
    }

    Info << "zoneSurfaceArea_ = " << zoneSurfaceArea_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcMassFluxSurface::read()
{
    IOdictionary dict
    (
        IOobject
        (
            "massFluxSurface_" + fieldName_ + "_" + faceZoneName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("molsZone", molsZone_);
    dict.readIfPresent("massZone", massZone_);

    dict.readIfPresent("averagingCounter", averagingCounter_);
}


void Foam::dsmcMassFluxSurface::write()
{
    if (mesh_.time().writeTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "massFluxSurface_"+fieldName_+"_"+faceZoneName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("molsZone", molsZone_);
        dict.add("massZone", massZone_);

        dict.add("averagingCounter", averagingCounter_);

        dict.regIOobject::writeObject
        (
            IOstreamOption(mesh_.time().writeFormat()),
            true
        );
    }
}


void Foam::dsmcMassFluxSurface::createField()
{}


void Foam::dsmcMassFluxSurface::calculateField()
{
    // Call this function every time-step before the state and flux objects are
    // cleaned
    ++sampleCounter_;

    if (sampleInterval_ <= sampleCounter_)
    {
        ++averagingCounter_;

        const List<scalarField>& molIdFlux = cloud_.tracker().parcelIdFlux();
        const List<scalarField>& massIdFlux = cloud_.tracker().massIdFlux();

        scalar molFlux = 0.0;
        scalar massFlux = 0.0;

        const faceZoneMesh& faceZones = mesh_.faceZones();
        const labelList& faces = faceZones[regionId_];

        for (const label faceI : faces)
        {
            vector nF = mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

            forAll(molIdFlux, id)
            {
                if (typeIds_.find(id) != -1)
                {
                    scalar RWF =
                        cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

                    molFlux +=
                        (molIdFlux[id][faceI]*cloud_.nParticle()*RWF*nF)
                      & fluxDirection_;
                    massFlux +=
                        (massIdFlux[id][faceI]*cloud_.nParticle()*RWF*nF)
                      & fluxDirection_;
                }
            }
        }

        molsZone_ += molFlux;
        massZone_ += massFlux;

        sampleCounter_ = 0;
    }

    const Time& runTime = mesh_.time();

    // Average measurement and calculate properties
    if (runTime.writeTime())
    {
        scalar molsZone = molsZone_;
        scalar massZone = massZone_;

        if (Pstream::parRun())
        {
            reduce(molsZone, sumOp<scalar>());
            reduce(massZone, sumOp<scalar>());
        }

        scalar averagingTime = averagingCounter_*timeVel_.mdTimeInterval().deltaT();

        molFluxZone_[timeIndex_] = molsZone/(averagingTime*zoneSurfaceArea_);
        massFluxZone_[timeIndex_] = massZone/(averagingTime*zoneSurfaceArea_);
        massFlowZone_[timeIndex_] = massZone/averagingTime;

        if (timeVel_.resetFieldsAtOutput())
        {
            molsZone_ = 0.0;
            massZone_ = 0.0;
            averagingCounter_ = 0.0;
        }

        if (averagingAcrossManyRuns_)
        {
            write();
        }

        ++timeIndex_;
    }
}


void Foam::dsmcMassFluxSurface::writeField()
{
    const Time& runTime = mesh_.time();

    if (runTime.writeTime())
    {
        timeIndex_ = 0;

        if (Pstream::master())
        {
            scalarField timeField(1);

            timeField[0] = mesh_.time().timeOutputValue();

            writeTimeData
            (
                casePath_,
                "faceFlux_" + faceZoneName_ + "_"+fieldName_ + "_N.xy",
                timeField,
                molFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceFlux_" + faceZoneName_ + "_" + fieldName_ + "_M.xy",
                timeField,
                massFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceMassFlowRate_" + faceZoneName_ + "_"+fieldName_ + ".xy",
                timeField,
                massFlowZone_,
                true
            );
        }
    }
}


void Foam::dsmcMassFluxSurface::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    dsmcField::updateProperties(dict);
}


// ************************************************************************* //
