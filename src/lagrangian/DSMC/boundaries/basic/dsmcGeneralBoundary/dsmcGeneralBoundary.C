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

#include "dsmcGeneralBoundary.H"
#include "fvMesh.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dsmcGeneralBoundary, 0);
defineRunTimeSelectionTable(dsmcGeneralBoundary, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcGeneralBoundary::dsmcGeneralBoundary
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcBoundaryBase(mesh, cloud, dict, "general"),
    faces_(),
    patchSurfaceArea_(0.0),
    cells_(),
    density_(0.0),
    velocity_(Zero),
    translationalTemperature_(0.0),
    rotationalTemperature_(0.0),
    vibrationalTemperature_(0.0),
    electronicTemperature_(0.0),
    accumulatedParcelsToInsert_()
{
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    // Initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    // Loop through all faces and set the boundary cells
    // - no conflict with parallelisation because the faces are unique

    for (label i = 0; i < patch.size(); ++i)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if (Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dsmcGeneralBoundary> Foam::dsmcGeneralBoundary::New
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("boundaryModel"));

    Info<< "Selecting dsmcGeneralBoundaryModel " << modelType << nl;

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

    return autoPtr<dsmcGeneralBoundary>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcGeneralBoundary::updateTime()
{}


void Foam::dsmcGeneralBoundary::computeParcelsToInsert
(
    const scalar& transT,
    const vector& velocity,
    const scalarList& numDen
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        const scalar& mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity & -sF/fA)/mostProbableSpeed;

            // From Bird eqn 4.22

            if (cloud_.axisymmetric())
            {
                scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

                accumulatedParcelsToInsert_[i][f] +=
                (
                    fA*numDen[i]*deltaT*mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
               /(2.0*sqrtPi*cloud_.nParticle()*RWF);
            }
            else
            {
                accumulatedParcelsToInsert_[i][f] +=
                (
                    fA*numDen[i]*deltaT*mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticle());
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::computeParcelsToInsert
(
    const scalarField& transT,
    const vectorField& velocity,
    const List<scalarField>& numDen
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    const scalar sqrtPi = sqrt(pi);

    // Compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        const scalar& mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT[f],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)
                /mostProbableSpeed;

            // From Bird eqn 4.22
            if (cloud_.axisymmetric())
            {
                scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

                accumulatedParcelsToInsert_[i][f] +=
                (
                    fA*numDen[i][f]*deltaT*mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticle()*RWF);
            }
            else
            {
                accumulatedParcelsToInsert_[i][f] +=
                (
                    fA*numDen[i][f]*deltaT*mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticle());
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::computeParcelsToInsert
(
    const scalar& transT,
    const vectorField& velocity,
    const scalar& numDen,
    const scalarField& molFractions
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar& fA = mag(sF);

            const scalar& pMass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    pMass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)
                                    /mostProbableSpeed;

            // From Bird eqn 4.22

            if (cloud_.axisymmetric())
            {
                const point& fC = cloud_.mesh().faceCentres()[faceI];
                scalar radius = fC.y();

                scalar RWF = 1.0;

                RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());

                accumulatedParcelsToInsert_[iD][f] +=
                    molFractions[iD]
                   *(
                        fA*numDen*deltaT*mostProbableSpeed
                       *(
                            exp(-sqr(sCosTheta))
                          + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                        )
                    )
                   /(2.0*sqrtPi*cloud_.nParticle()*RWF);
            }
            else
            {
                accumulatedParcelsToInsert_[iD][f] +=
                    molFractions[iD]*
                    (
                        fA*numDen*deltaT*mostProbableSpeed
                       *(
                            exp(-sqr(sCosTheta))
                          + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                        )
                    )
                   /(2.0*sqrtPi*cloud_.nParticle());
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::computeParcelsToInsert
(
    const scalarField& transT,
    const vectorField& velocity,
    const scalarField& numDen,
    const scalarField& molFractions
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            const scalar fA = mag(sF);

            scalar mass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT[f],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)/mostProbableSpeed;

            // From Bird eqn 4.22
//             scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faces_[f]]);
// 
//             accumulatedParcelsToInsert_[iD][f] +=
//                 molFractions[iD]
//                *(
//                     fA*numDen[f]*deltaT*mostProbableSpeed
//                    *(
//                         exp(-sqr(sCosTheta))
//                       + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
//                     )
//                 )
//                /(2.0*sqrtPi*cloud_.nParticle()*RWF);
            
            if (cloud_.axisymmetric())
            {
                const point& fC = cloud_.mesh().faceCentres()[faces_[f]];
                scalar radius = fC.y();

                scalar RWF = 1.0;

                RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());

                accumulatedParcelsToInsert_[iD][f] +=
                    molFractions[iD]
                   *(
                        fA*numDen[f]*deltaT*mostProbableSpeed
                       *(
                            exp(-sqr(sCosTheta))
                          + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                        )
                    )
                   /(2.0*sqrtPi*cloud_.nParticle()*RWF);
            }
            else
            {
                accumulatedParcelsToInsert_[iD][f] +=
                    molFractions[iD]*
                    (
                        fA*numDen[f]*deltaT*mostProbableSpeed
                       *(
                            exp(-sqr(sCosTheta))
                          + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                        )
                    )
                   /(2.0*sqrtPi*cloud_.nParticle());
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::computeParcelsToInsert
(
    const scalar& transT,
    const vectorField& velocity,
    const List<scalarField>& numDen
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            const scalar fA = mag(sF);

            scalar mass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA )/mostProbableSpeed;

            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faces_[f]]);

            accumulatedParcelsToInsert_[iD][f] +=
            (
                fA*numDen[iD][f]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
           /(2.0*sqrtPi*cloud_.nParticle()*RWF);
        }
    } 
}


void Foam::dsmcGeneralBoundary::insertParcels
(
    const scalar& transT,
    const scalar& rotT,
    const scalar& vibT,
    const scalar& elecT,
    const vector& velocity
)
{
    Random& rndGen = cloud_.rndGen();

    labelField parcelsInserted(typeIds_.size(), Zero);
    labelField parcelsToAdd(typeIds_.size(), Zero);

    // insert pacels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets =
            polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);

        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non - flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        forAll(typeIds_, m)
        {
            const label typeId = typeIds_[m];

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

            const scalar mass = cloud_.constProps(typeId).mass();

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

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        transT,
                        mass
                    )
                );

                scalar sCosTheta = (velocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5
                   *(
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if (mag(velocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P =
                                2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()
                   *transT/mass)
                   *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                  + (t1 & velocity)*t1
                  + (t2 & velocity)*t2
                  + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    rotT,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    vibT,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    elecT,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );

                label newParcel = 1;

                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += 1e-5*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

                parcelsInserted[m] += 1.0;
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::insertParcels
(
    const scalarField& transT,
    const scalarField& rotT,
    const scalarField& vibT,
    const scalarField& elecT,
    const vectorField& velocity
)
{
    Random& rndGen = cloud_.rndGen();

    labelField parcelsInserted(typeIds_.size(), 0);
    labelField parcelsToAdd(typeIds_.size(), 0);

    // Insert parcels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);
        const vector& faceVelocity = velocity[f];
        const scalar faceTranslationalTemperature = transT[f];
        const scalar faceRotationalTemperature = rotT[f];
        const scalar faceVibrationalTemperature = vibT[f];
        const scalar faceElectronicTemperature = elecT[f];

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);

        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non - flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        forAll(typeIds_, m)
        {
            const label typeId = typeIds_[m];

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

            scalar mass = cloud_.constProps(typeId).mass();

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

                // Velocity generation

                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTranslationalTemperature,
                        mass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5
                   *(
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                // Select a velocity using Bird eqn 12.5
                do
                {
                    uNormalThermal =
                        randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                    uNormal = uNormalThermal + sCosTheta;

                    if (uNormal < 0.0)
                    {
                        P = -1;
                    }
                    else
                    {
                        P =
                            2.0*uNormal/uNormProbCoeffA
                           *exp(uNormProbCoeffB - sqr(uNormalThermal));
                    }

                } while (P < rndGen.sample01<scalar>());

                vector U =
                    sqrt(physicoChemical::k.value()
                       *faceTranslationalTemperature/mass)
                       *(
                            rndGen.GaussNormal<scalar>()*t1
                          + rndGen.GaussNormal<scalar>()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceRotationalTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceVibrationalTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceElectronicTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );


                label newParcel = 1;

                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += 1e-5*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

                parcelsInserted[m] += 1.0;
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::insertParcels
(
    const scalar& transT,
    const vectorField& velocity
)
{
    Random& rndGen = cloud_.rndGen();

    label nTotalParcelsAdded = 0;
    label nTotalParcelsToBeAdded = 0;

    // Loop over all species
    forAll(accumulatedParcelsToInsert_, iD)
    {
        // Loop over all faces of the patch
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& faceVelocity = velocity[f];
            const scalar faceTemperature = transT;
            const label faceI = faces_[f];
            const label cellI = cells_[f];
            const vector& fC = mesh_.faceCentres()[faceI];
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            scalar fA = mag(sF);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh_).mag()/fA
                    + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non - flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //  Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            /* -------------------------------------------------------------*/

            // generate Poisson distributed random number of particles to insert
            // see Tysanner & Garcia, Int. J. Numer. Meth. Fluids 2050; 00:1-12
            // this eliminates non - equilibrium behaviour that does not exist
            // in the corresponding physical system

            label k = 0;

            const scalar target = exp(-accumulatedParcelsToInsert_[iD][f]);

            scalar p = rndGen.sample01<scalar>();

            while (p > target)
            {
                p *= rndGen.sample01<scalar>();
                k += 1;
            }

            label nParcelsToInsert = k;

            /* ----------------------------------------------------------*/
            
            accumulatedParcelsToInsert_[iD][f] -= nParcelsToInsert;
            

            nTotalParcelsToBeAdded += nParcelsToInsert;

            const label typeId = typeIds_[iD];
            const scalar pMass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; i++)
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

                point pt(faceTetIs.faceTri(mesh_).randomPoint(rndGen));

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTemperature,
                        pMass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5*
                    (
                        1.0
                        + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if (mag(faceVelocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                                *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()*faceTemperature/pMass)
                   *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                  + (t1 & faceVelocity)*t1
                  + (t2 & faceVelocity)*t2
                  + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );
                
                label newParcel = 1;

                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                pt += (n*SMALL);

                cloud_.addNewParcel
                (
                    pt,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

                ++nTotalParcelsAdded;
            }
        }
    }
}

void Foam::dsmcGeneralBoundary::insertParcels
(
    const scalarField& transT,
    const vectorField& velocity
)
{
    Random& rndGen = cloud_.rndGen();

    label nTotalParcelsAdded = 0;
    label nTotalParcelsToBeAdded = 0;

    labelField parcelsInserted(typeIds_.size(), 0);

    // loop over all species
    forAll(accumulatedParcelsToInsert_, iD)
    {
        // loop over all faces of the patch
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& faceVelocity = velocity[f];
            const scalar faceTemperature = transT[f];
            const label faceI = faces_[f];
            const label cellI = cells_[f];
            const vector& fC = mesh_.faceCentres()[faceI];
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            scalar fA = mag(sF);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

            //Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh_).mag()/fA
                    + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }


            //Force the last area fraction value to 1.0 to avoid any
            //rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            //Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //Wall tangential unit vector. Use the direction between the
            //face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            //Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            label nParcelsToInsert = label(accumulatedParcelsToInsert_[iD][f]);

            if ((nParcelsToInsert - accumulatedParcelsToInsert_[iD][f]) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            // Note: remainder has been set
            accumulatedParcelsToInsert_[iD][f] -= nParcelsToInsert;

            nTotalParcelsToBeAdded += nParcelsToInsert;

            const label typeId = typeIds_[iD];
            const scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; ++i)
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

                point p (faceTetIs.faceTri(mesh_).randomPoint(rndGen));

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTemperature,
                        mass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5
                   *(
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if (mag(faceVelocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()*faceTemperature/mass)
                   *(
                        rndGen.GaussNormal<scalar>()*t1
                      + rndGen.GaussNormal<scalar>()*t2
                    )
                  + (t1 & faceVelocity)*t1
                  + (t2 & faceVelocity)*t2
                  + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label newParcel = 1;

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );

                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += 1e-5*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

                ++nTotalParcelsAdded;
                ++parcelsInserted[iD];
            }
        }
    }
    
}


const Foam::labelList& Foam::dsmcGeneralBoundary::controlPatch() const
{
    return faces_;
}


const Foam::labelList& Foam::dsmcGeneralBoundary::controlZone() const
{
    return cells_;
}


Foam::scalar Foam::dsmcGeneralBoundary::density() const
{
    return density_;
}


Foam::scalar& Foam::dsmcGeneralBoundary::density()
{
    return density_;
}


const Foam::vector& Foam::dsmcGeneralBoundary::velocity() const
{
    return velocity_;
}


Foam::vector& Foam::dsmcGeneralBoundary::velocity()
{
    return velocity_;
}


Foam::scalar Foam::dsmcGeneralBoundary::translationalTemperature() const
{
    return translationalTemperature_;
}


Foam::scalar& Foam::dsmcGeneralBoundary::translationalTemperature()
{
    return translationalTemperature_;
}


Foam::scalar Foam::dsmcGeneralBoundary::rotationalTemperature() const
{
    return rotationalTemperature_;
}


Foam::scalar& Foam::dsmcGeneralBoundary::rotationalTemperature()
{
    return rotationalTemperature_;
}

Foam::scalar Foam::dsmcGeneralBoundary::vibrationalTemperature() const
{
    return vibrationalTemperature_;
}


Foam::scalar& Foam::dsmcGeneralBoundary::vibrationalTemperature()
{
    return vibrationalTemperature_;
}


Foam::scalar Foam::dsmcGeneralBoundary::electronicTemperature() const
{
    return electronicTemperature_;
}


Foam::scalar& Foam::dsmcGeneralBoundary::electronicTemperature()
{
    return electronicTemperature_;
}


// ************************************************************************* //
