/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

Measures overall temperature, including vibrational temperature, for a single 
species gas or a gas mixture and writes the results to a volume scalar field 
that can be viewed in Paraview.

Translational, rotatational and vibrational temperature field will also be 
written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements
class and are also written.

\*---------------------------------------------------------------------------*/

#include "solidVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidVolFields, 0);

addToRunTimeSelectionTable(solidField, solidVolFields, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidVolFields::solidVolFields
(
    Time& t,
    const polyMesh& mesh,
    solidParcelCloud& cloud,
    const dictionary& dict
)
:
    solidField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    n_(),
    t1_(),
    t2_(),
    sampleInterval_(1),
    sampleCounter_(0),
    nTimeSteps_(0),
    fieldName_(propsDict_.lookup("fieldName")),
    typeIds_(),
    nParcelMean_
    (
        IOobject
        (
            "nParcelMean_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    ),
    q_
    (
        IOobject
        (
            "q_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    UMean_
    (
        IOobject
        (
            "UMean_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    nParcelsCum_(mesh_.nCells(), 0.0),
    nCum_(mesh_.nCells(), 0.0),
    mCum_(mesh_.nCells(), 0.0),
    momentumCum_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    rhoNBF_(),
    rhoMBF_(),
    qBF_(),
    momentumBF_(),
    fDBF_(),
    averagingAcrossManyRuns_(false)
{

    // standard to reading typeIds ------------ 
    typeIds_ = cloud_.getTypeIDs(propsDict_);
    
    boundaryCells_.setSize(mesh.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];

        boundaryCells_[p].setSize(patch.size());

        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }
    
    // initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());

    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());
        
    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        
        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        fDBF_[j].setSize(patch.size(), vector::zero);
        
        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }
    
    calculateWallUnitVectors();
    
    sampleInterval_ = propsDict_.lookupOrDefault("sampleInterval", 1);

    averagingAcrossManyRuns_ = 
        propsDict_.lookupOrDefault<bool>("averagingAcrossManyRuns", false);
    
    if(averagingAcrossManyRuns_)
    {
        // read in stored data from dictionary
        if (averagingAcrossManyRuns_)
        {
            Info<< "averagingAcrossManyRuns initiated." << nl << endl;
            readIn();
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidVolFields::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "volFields_"+fieldName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );
    
    dict.readIfPresent("nTimeSteps", nTimeSteps_);

    // cumulative values
    dict.readIfPresent("nParcelsCum", nParcelsCum_);
    dict.readIfPresent("nCum", nCum_);
    dict.readIfPresent("mCum", mCum_);
    dict.readIfPresent("momentumCum", momentumCum_);

    dict.readIfPresent("rhoNBF", rhoNBF_);
    dict.readIfPresent("rhoMBF", rhoMBF_);
    dict.readIfPresent("qBF", qBF_); 
    dict.readIfPresent("momentumBF", momentumBF_);
    dict.readIfPresent("fDBF", fDBF_); 
}

void solidVolFields::writeOut()
{
    if (mesh_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFields_"+fieldName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("nTimeSteps", nTimeSteps_);
        // cumulative values
        dict.add("nCum", nCum_);
        dict.add("mCum", mCum_);
        dict.add("nParcelsCum", nParcelsCum_);
        dict.add("momentumCum", momentumCum_);

        dict.add("rhoNBF", rhoNBF_);
        dict.add("rhoMBF", rhoMBF_);
        dict.add("qBF", qBF_); 
        dict.add("momentumBF", momentumBF_);
        dict.add("fDBF", fDBF_); 
        
        IOstream::streamFormat fmt = mesh_.time().writeFormat();
        IOstream::versionNumber ver = mesh_.time().writeVersion();
        IOstream::compressionType cmp = mesh_.time().writeCompression();
    
        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

void solidVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];
       
        if(isA<wallPolyPatch>(pPatch))
        {           
            const vectorField& fC = pPatch.faceCentres();
           
            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei]/mag(pPatch.faceAreas()[facei]);
               
                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei] - 
                            mesh_.points()[mesh_.faces()[pPatch.start() 
                            + facei][0]];
                t1_[patchi][facei] /= mag(t1_[patchi][facei]);
               
                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei] /= mag(t2_[patchi][facei]);
            }
        }
    }
}

//- initial conditions
void solidVolFields::createField()
{
    Info << "Initialising solidVolFields field" << endl;
}


void solidVolFields::calculateField()
{
    ++sampleCounter_;
    
    if(sampleInterval_ <= sampleCounter_)
    {
        ++nTimeSteps_;
        
        forAllConstIters(cloud_, iter)
        {
            const solidParcel& p = iter();
            const label iD = findIndex(typeIds_, p.typeId());

            if(iD != -1)
            {
                const label cell = p.cell();
                const scalar nParticles = cloud_.nParticle();
                const scalar mass = p.mass();     
                const vector& U = p.U();
                                
                nParcelsCum_[cell] += 1.0;

                if(cloud_.axisymmetric())
                {
                    const point& cC = cloud_.mesh().cellCentres()[cell];
                    
                    scalar radius = cC.y();
                    
                    scalar RWF = 1.0;
                    
                    RWF = 1.0 + 
                        cloud_.maxRWF()*(radius/cloud_.radialExtent());
                    nCum_[cell] += RWF*nParticles;
                    mCum_[cell] += mass*RWF*nParticles;
                
                    momentumCum_[cell] += (mass*(U)*RWF*nParticles);
                    //rhoNMeanXnParticle_[cell] += (RWF*nParticle);
                    //rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                }
                else
                {
                    nCum_[cell] += nParticles;
                    mCum_[cell] += mass*nParticles;
                
                    momentumCum_[cell] += (mass*(U)*nParticles);
                }
                //nCum_[cell] += nParticles;
                //mCum_[cell] += mass*nParticles;
                
                //momentumCum_[cell] += (mass*(U)*nParticles);
            }
        }
        
        
        // Obtain collision quality measurements
        
        //- Obtain boundary measurements
        forAll(typeIds_, i)
        {
            const label typeId = typeIds_[i];
            
            forAll(mesh_.boundaryMesh(), j)
            {                
                forAll(mesh_.boundaryMesh()[j], k)
                {
                    rhoNBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[typeId][j][k];
                    rhoMBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoMBF()[typeId][j][k];
                    momentumBF_[j][k] += cloud_.boundaryFluxMeasurements().momentumBF()[typeId][j][k];
                    fDBF_[j][k] += cloud_.boundaryFluxMeasurements().fDBF()[typeId][j][k];
                    qBF_[j][k] += cloud_.boundaryFluxMeasurements().qBF()[typeId][j][k];
                }
            }
        }

        sampleCounter_ = 0;
    }
    
    if(mesh_.time().writeTime())//- output of writeTime() is writeTime_ 
    {
        const scalar nAvTimeSteps = nTimeSteps_;

            forAll(nParcelsCum_, cell)
            {      
                if(nParcelsCum_[cell] > SMALL)
                {  
                    const scalar cellVolume = mesh_.cellVolumes()[cell];
                    nParcelMean_[cell] = nParcelsCum_[cell]/nAvTimeSteps;
                    rhoN_[cell] = (nCum_[cell])/(nAvTimeSteps*cellVolume);
                    rhoM_[cell] = (mCum_[cell])/(nAvTimeSteps*cellVolume);
                    UMean_[cell] = momentumCum_[cell] / (mCum_[cell]);

                }
                else
                {
                    // not zero so that weighted decomposition still works
                    nParcelMean_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                } 
            }
            
            // computing boundary measurements
            forAll(rhoNBF_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                const bool isWall = isA<wallPolyPatch>(patch);
                const bool isNonEmptyNonCyclic = isA<polyPatch>(patch)
                    && !isA<emptyPolyPatch>(patch)
                    && !isA<cyclicPolyPatch>(patch);

                if(isWall)
                {                               
                    forAll(patch, k)
                    {                 
                        const label celli = boundaryCells_[j][k];    
                        // Note: do not use the nParticles value that includes
                        // the RWF here. This is wrong because boundary
                        // measurements are performend during move steps. Hence
                        // radial weighting (if simulation is axisymmetric or
                        // spherical) is performed after the measurement. That
                        // is why the parcel RWFs are already included during
                        // the measurement step and we only need the FNUM value
                        // here.
                        //const scalar nParticles = cloud_.coordSystem().dtModel().nParticles(j, k);
                        const scalar nParticles = cloud_.nParticle();

                        nParcelMean_.boundaryFieldRef()[j][k] = nParcelMean_[celli];
                        rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*nParticles/nAvTimeSteps;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*nParticles/nAvTimeSteps;
                        
                        if(rhoM_.boundaryFieldRef()[j][k] > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]/rhoMBF_[j][k];
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                        }

                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/(nAvTimeSteps);
                    }


                }
                else if (isNonEmptyNonCyclic)
                {
                    //- Loop over all boundary faces and set zeroGradient conditions
                    forAll(boundaryCells_[j], k)
                    {
                        const label celli = boundaryCells_[j][k];
                        
                        nParcelMean_.boundaryFieldRef()[j][k] = nParcelMean_[celli];
                        rhoN_.boundaryFieldRef()[j][k] = rhoN_[celli];
                        rhoM_.boundaryFieldRef()[j][k] = rhoM_[celli];
                        
                        UMean_.boundaryFieldRef()[j][k] = UMean_[celli];
                    }
                }
            }
            
            UMean_.write();
            q_.write();
        
        //- reset
        if(timeDM_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0;
            
            forAll(nParcelsCum_, c)
            {
                nParcelsCum_[c] = 0.0;
                
                nCum_[c] = 0.0;
                mCum_[c] = 0.0;

                momentumCum_[c] = vector::zero;
            }
            
            // reset boundary information
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                qBF_[j] = 0.0;
                momentumBF_[j] = vector::zero;
                fDBF_[j] = vector::zero;
            }
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
        
    }
}


//- write field
void solidVolFields::writeField()
{}

volScalarField& solidVolFields::surfaceHeatFlux()
{
    return q_;
}

const word& solidVolFields::fieldName() const
{
    return fieldName_;
}

void solidVolFields::updateProperties(const dictionary& dict)
{
    //- the main properties should be updated first
    solidField::updateProperties(dict);
    
    propsDict_ = dict.subDict(typeName + "Properties");
    
    Info << "Update Field Properties !"<< endl;
    
    propsDict_.readIfPresent
    (
        "averagingAcrossManyRuns",
        averagingAcrossManyRuns_
    );
    if (averagingAcrossManyRuns_)
    {
        Info << "averagingAcrossManyRuns initiated." << endl;
    }

}


} // End namespace Foam

// ************************************************************************* //

