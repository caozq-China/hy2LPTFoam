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

#include "solidGeneralBoundary.H"
#include "fvMesh.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "solidParticleCouplingCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidGeneralBoundary, 0);

defineRunTimeSelectionTable(solidGeneralBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidGeneralBoundary::solidGeneralBoundary
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidBoundaryBase(mesh, spc, dict, "general"),
    faces_(),
    patchSurfaceArea_(0.0),
    cells_(),
//     phaseStates_(),
//     CzRatios_(),
    accumulatedParcelsToInsert_()
{
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    for(label i = 0; i < patch.size(); ++i)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidGeneralBoundary> solidGeneralBoundary::New
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
{
    word modelType
    (
        dict.get<word>("boundaryModel")
    );

    Info<< "Selecting solidGeneralBoundaryModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "boundaryModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<solidGeneralBoundary>
    (
        cstrIter()(mesh, spc, dict)
    );
}

void solidGeneralBoundary::updateTime()
{}

void Foam::solidGeneralBoundary::computeParcelsToInsert
(
    const vector& velocity,
    const scalarList& numDen
)
{
    const scalar deltaT = mesh_.time().deltaTValue();

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);
            
            const vector sFnormal = sF/mag(sF);

            // From Bird eqn 4.22
            
            if(spc_.dsmcCloudReference()->axisymmetric())
            {
                scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().faceCentres()[faceI]);
                
                accumulatedParcelsToInsert_[i][f] += 
                    fA*numDen[i]*deltaT*(velocity & (-sFnormal))/(spc_.nSolidParticles()*RWF);
            }
            else
            {
                accumulatedParcelsToInsert_[i][f] += 
                    fA*numDen[i]*deltaT*(velocity & (-sFnormal))/spc_.nSolidParticles();
            }
        }
    }
}

void Foam::solidGeneralBoundary::insertParcels
(
    const scalar& temperature,
    const vector& velocity,
//     const scalarList& phaseStates,
    const scalarList& CzRatios
)
{
    
    Random& rndGen = spc_.rndGenS();
    
    labelField parcelsInserted(typeIDs_.size(), Zero);
    labelField parcelsToAdd(typeIDs_.size(), Zero);
    // insert solid particles
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
//         const vector fC = mesh_.faceCentres()[faceI];
        const vector sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        //         Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);
        
        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        //         Force the last area fraction value to 1.0 to avoid any
        //         rounding/non-flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        forAll(typeIDs_, m)
        {
            const label typeIdSolid = typeIDs_[m];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[m][f];
            
            // Number of whole particles to insert
            label nI = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
            {
                ++nI;
            }
            
            faceAccumulator -= nI;
            parcelsToAdd[m] += nI;
//             solidToAdd[m] += nI;

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

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);
                
                label newSolidParticle = 1;
                
                scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().cellCentres()[cellI]);
                
                scalar particleDiameter = 0.0;
                
                label phaseState = 0;
                
                if(spc_.solidPhaseChange().active())
                {
                    //- particle phase check
                    if(temperature< spc_.constSolidProps(m).Tf() && CzRatios[m] ==0.0)
                    {
                        //- pure solid phase
                        phaseState = 0;
                    }
                    else if((temperature < spc_.constSolidProps(m).Tf() && CzRatios[m] ==1.0) || (temperature<spc_.constSolidProps(m).Tm() && (CzRatios[m]<1.0 && CzRatios[m]>0.0)))
                    {
                        //- the core is liquid
                        phaseState = 1;
                    }
                    else if(temperature >= spc_.constSolidProps(m).Tm() && CzRatios[m] < 1.0 && CzRatios[m] > 0.0 )
                    {
                        //- the core is solid
                        phaseState = 2;
                    }
                    else if(temperature >= spc_.constSolidProps(m).Tm() && CzRatios[m] ==0.0)
                    {
                        //- pure liquid phase
                        phaseState = 3;
                    }
                    
                    
                    //- particle diameter correction
                    if(phaseState == 0 || phaseState == 1)
                    {
                        //- particle diameter will not be corrected through density correction equation when 
                        //- particle temperature is lower than melting temperature.
                        particleDiameter = spc_.constSolidProps(m).d();
                    }
                    else if(phaseState == 2)
                    {
                        scalar rhoLiquid = spc_.particleLiquidDensityCorrection
                                                (
                                                    spc_.materialList()[m],
                                                    temperature
                                                );
                        
                        particleDiameter = pow(
                                        spc_.constSolidProps(m).massSphere()/
                                        (pi*(pow(CzRatios[m],3.0)*spc_.constSolidProps(m).rho()+(1.0-pow(CzRatios[m],3.0))*rhoLiquid)/6.0)
                                        ,1.0/3.0
                                    );
                    }
                    else if(phaseState == 3)
                    {
                        scalar rhoLiquid = spc_.particleLiquidDensityCorrection
                                            (
                                                spc_.materialList()[m],
                                                temperature
                                            );
                    
                        particleDiameter = pow((6.0*spc_.constSolidProps(m).massSphere()/(pi*rhoLiquid)),1.0/3.0);
                    }
                }
                else
                {
                    particleDiameter = spc_.constSolidProps(m).d();
                    
                }
                
                //- just for MPPIC validation
//                 scalar PHI = 360*rndGen.scalar01();
//                 USolid_.y() = sin(PHI);
//                 USolid_.z() = cos(PHI);
                
                
                spc_.addNewSolidParticle
                (
                    p,
                    velocity,
                    (vector::zero),
                    (vector::zero),
                    (vector::zero),
                    particleDiameter,
                    temperature,
                    RWF,
                    CzRatios[m],
                    cellI,
                    -1,
                    typeIdSolid,
                    phaseState,
                    newSolidParticle
                );
                
                parcelsInserted[m] += 1.0;

            }
        }
    }
}

void Foam::solidGeneralBoundary::insertParcels
(
    const scalar& temperature,
    vector& velocity,
    const scalar& nozzleRadius,
    const scalar& maxOffAxisAngle,
//     const scalarList& phaseStates,
    const scalarList& CzRatios
)
{
    Random& rndGen = spc_.rndGenS();
    
    labelField parcelsInserted(typeIDs_.size(), Zero);
    labelField parcelsToAdd(typeIDs_.size(), Zero);
    // insert solid particles
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
//         const vector fC = mesh_.faceCentres()[faceI];
        const vector sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        //         Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);
        
        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        //         Force the last area fraction value to 1.0 to avoid any
        //         rounding/non-flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        forAll(typeIDs_, m)
        {
            const label typeIdSolid = typeIDs_[m];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[m][f];
            
            // Number of whole particles to insert
            label nI = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
            {
                ++nI;
            }
            
            faceAccumulator -= nI;
            parcelsToAdd[m] += nI;
//             solidToAdd[m] += nI;

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

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);
                
                //- radial direction: y
                //- axial direction: x
                //- calculate the off-axis angle and change its unit from degree to radian
                scalar offaxisAngle = (p.y()/nozzleRadius)*maxOffAxisAngle*pi/180;
                
                scalar magUxy = sqrt(velocity.x()*velocity.x()+velocity.y()*velocity.y());
                
                velocity.x() = magUxy*cos(offaxisAngle);
                velocity.y() = magUxy*sin(offaxisAngle);
                
                label newSolidParticle = 1;
                
                scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().cellCentres()[cellI]);
                
                scalar particleDiameter = 0.0;
                label phaseState = 0;
                
                if(spc_.solidPhaseChange().active())
                {
                    //- particle phase check
                    if(temperature< spc_.constSolidProps(m).Tf() && CzRatios[m] ==0.0)
                    {
                        //- pure solid phase
                        phaseState = 0;
                    }
                    else if((temperature < spc_.constSolidProps(m).Tf() && CzRatios[m] ==1.0) || (temperature<spc_.constSolidProps(m).Tm() && (CzRatios[m]<1.0 && CzRatios[m]>0.0)))
                    {
                        //- the core is liquid
                        phaseState = 1;
                    }
                    else if(temperature >= spc_.constSolidProps(m).Tm() && CzRatios[m] < 1.0 && CzRatios[m] > 0.0 )
                    {
                        //- the core is solid
                        phaseState = 2;
                    }
                    else if(temperature >= spc_.constSolidProps(m).Tm() && CzRatios[m] ==0.0)
                    {
                        //- pure liquid phase
                        phaseState = 3;
                    }
                    
                    
                    //- particle diameter correction
                    if(phaseState == 0 || phaseState == 1)
                    {
                        //- particle diameter will not be corrected through density correction equation when 
                        //- particle temperature is lower than melting temperature.
                        particleDiameter = spc_.constSolidProps(m).d();
                    }
                    else if(phaseState == 2)
                    {
                        scalar rhoLiquid = spc_.particleLiquidDensityCorrection
                                                (
                                                    spc_.materialList()[m],
                                                    temperature
                                                );
                        
                        particleDiameter = pow(
                                        spc_.constSolidProps(m).massSphere()/
                                        (pi*(pow(CzRatios[m],3.0)*spc_.constSolidProps(m).rho()+(1.0-pow(CzRatios[m],3.0))*rhoLiquid)/6.0)
                                        ,1.0/3.0
                                    );
                    }
                    else if(phaseState == 3)
                    {
                        scalar rhoLiquid = spc_.particleLiquidDensityCorrection
                                            (
                                                spc_.materialList()[m],
                                                temperature
                                            );
                    
                        particleDiameter = pow((6.0*spc_.constSolidProps(m).massSphere()/(pi*rhoLiquid)),1.0/3.0);
                    }
                }
                else
                {
                    particleDiameter = spc_.constSolidProps(m).d();
                    
                }
                
                //- just for MPPIC validation
//                 scalar PHI = 360*rndGen.scalar01();
//                 USolid_.y() = sin(PHI);
//                 USolid_.z() = cos(PHI);
                
                
                spc_.addNewSolidParticle
                (
                    p,
                    velocity,
                    (vector::zero),
                    (vector::zero),
                    (vector::zero),
                    particleDiameter,
                    temperature,
                    RWF,
                    CzRatios[m],
                    cellI,
                    -1,
                    typeIdSolid,
                    phaseState,
                    newSolidParticle
                );
                
                parcelsInserted[m] += 1.0;

            }
        }
    }
}


const labelList& solidGeneralBoundary::controlPatch() const
{
    return faces_;
}

const labelList& solidGeneralBoundary::controlZone() const
{
    return cells_;
}

const vector& solidGeneralBoundary::velocity() const
{
    return velocity_;
}


vector& solidGeneralBoundary::velocity()
{
    return velocity_;
}

const scalar& solidGeneralBoundary::temperature() const
{
    return temperature_;
}


scalar& solidGeneralBoundary::temperature()
{
    return temperature_;
}


} // End namespace Foam

// ************************************************************************* //
