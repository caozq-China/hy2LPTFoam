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

\*---------------------------------------------------------------------------*/

#include "solidMassLoadingRatioInflowPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "solidParcelCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(solidMassLoadingRatioInflowPatch, 0);

addToRunTimeSelectionTable(solidGeneralBoundary, solidMassLoadingRatioInflowPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidMassLoadingRatioInflowPatch::solidMassLoadingRatioInflowPatch
(
    const polyMesh& mesh,
    solidParcelCloud& spc,
    const dictionary& dict
)
:
    solidGeneralBoundary(mesh, spc, dict),
    solidPropsDict_(dict.subDict(typeName + "Properties")),
    massLoadingRatio_(solidPropsDict_.lookup("massLoadingRatio")),
    parcelMass_(massLoadingRatio_.size(),0.0),
    numberDensities_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    // setProperties();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void solidMassLoadingRatioInflowPatch::initialConfiguration()
{
    setProperties();
}

void solidMassLoadingRatioInflowPatch::calculateProperties()
{}

void solidMassLoadingRatioInflowPatch::controlParcelsBeforeMove()
{
    const scalar deltaT = mesh_.time().deltaTValue();

    const vectorField& UcB = spc_.UFilter().boundaryField()[patchId_];
    const scalarField& TcB = spc_.TtrFilter().boundaryField()[patchId_];

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);
            const vector sFnormal = sF/mag(sF);

            //accumulatedParcelsToInsert_[i][f] +=
            //        fA*numberDensities_[i][f]*deltaT*(UcB[f] & (-sFnormal))/spc_.nParticles(patchId_,f);
            if(spc_.axisymmetric())
            {
                const point fC = spc_.mesh().faceCentres()[faceI];
                scalar radius = fC.y();
                
                scalar RWF = 1.0;
                    
                RWF = 1.0 + spc_.maxRWF()*(radius/spc_.radialExtent());
                accumulatedParcelsToInsert_[i][f] +=
                    fA*numberDensities_[i][f]*deltaT*(UcB[f] & (-sFnormal))/(spc_.nParticle()*RWF);
            }
            else
            {
                accumulatedParcelsToInsert_[i][f] +=
                    fA*numberDensities_[i][f]*deltaT*(UcB[f] & (-sFnormal))/spc_.nParticle();
            }
        }
    }

    Random& rndGen = spc_.rndGen();
    
    // insert solid particles
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), 0.0);

        scalar previousCumulativeSum = 0.0;

        forAll(faceTets, triI)
        {
            const tetIndices faceTetIs = faceTets[triI];

            cTriAFracs[triI] =
                faceTetIs.faceTri(mesh_).mag()/fA
                + previousCumulativeSum;

            previousCumulativeSum = cTriAFracs[triI];
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non-flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;


        forAll(typeIDs_, m)
        {
            const label typeId = typeIDs_[m];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[m][f];

            // Number of whole particles to insert
            label nI = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
            {
                nI++;
            }

            faceAccumulator -= nI;

            for (label i = 0; i < nI; ++i)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point position = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                //const scalar RWF = spc_.coordSystem().RWF(cellI);
                scalar RWF = 1.0;
                    
                if(spc_.axisymmetric())
                {                      
                    const point cC = spc_.mesh().cellCentres()[cellI];
                    scalar radius = cC.y();

                    RWF = 1.0 + spc_.maxRWF()*(radius/spc_.radialExtent());
                }
                
                spc_.addNewParcel
                (
                    mesh_,
                    spc_.constProps(typeId),
                    position,
                    UcB[f],
                    TcB[f],
                    RWF,
                    cellI,
                    faces_[f],
                    faceTetIs.tetPt(),
                    typeId,
                    -1    
                );

            }
        }
    }
}

void solidMassLoadingRatioInflowPatch::controlParcelsBeforeCollisions()
{}

void solidMassLoadingRatioInflowPatch::controlParcelsAfterCollisions()
{}

void solidMassLoadingRatioInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void solidMassLoadingRatioInflowPatch::updateProperties(const dictionary& Dict)
{
    //- the main properties should be updated first
    solidGeneralBoundary::updateProperties(Dict);
}



void solidMassLoadingRatioInflowPatch::setProperties()
{
    //  read in the type ids
    typeIDs_ = spc_.getTypeIDs(solidPropsDict_);

    const scalarField& rhocB = spc_.rhoFilter().boundaryField()[patchId_];

    forAll(typeIDs_,i)
    {
        parcelMass_[i] = spc_.constProps(i).rho0()*pi*pow(spc_.constProps(i).d0(),3)/6;
    }

    //- number density initialization
    numberDensities_.clear();

    numberDensities_.setSize(typeIDs_.size());

    forAll(numberDensities_, i)
    {
        numberDensities_[i].setSize(faces_.size(),Zero);
        forAll(numberDensities_[i],f)
        {
            numberDensities_[i][f] = massLoadingRatio_[i]*rhocB[f]/parcelMass_[i];
        }
        
    }
    
    // set the accumulator
    accumulatedParcelsToInsert_.setSize(typeIDs_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(),Zero);
    }
}



} // End namespace Foam

// ************************************************************************* //
