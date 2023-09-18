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

#include "xParticleMeshFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(xParticleMeshFill, 0);

addToRunTimeSelectionTable(configuration, xParticleMeshFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
xParticleMeshFill::xParticleMeshFill
(
    xParticleCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    configuration(cloud, dict)
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void xParticleMeshFill::setInitialConfiguration()
{
   
    
    Info<< nl << "Initialising  extraterrestrial particles" << endl;

    // const scalar interphaseInitialRelativeSpeed = solidInitialiseDict_.get<scalar>("interphaseInitialRelativeSpeed");

    // const scalar interparticleInitialRelativeSpeed = solidInitialiseDict_.get<scalar>("interparticleInitialRelativeSpeed");
       

    const dictionary& numberDensitiesDict
    (
        initialiseDict_.subDict("numberDensities")
    );
    //- number density
    List<word> particles(numberDensitiesDict.toc());

    Field<scalar> numberDensities(particles.size());

    forAll(particles, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(particles[i]);
    }
    
    const dictionary& velocityDict
    (
        initialiseDict_.subDict("velocity")
    );
    Field<vector> velocities(particles.size(),vector::zero);
    
    forAll(particles, i)
    {
        velocities[i] = velocityDict.get<vector>(particles[i]);
    }

//     const dictionary& velocityDict
//     (
//         initialiseDict_.subDict("velocity")
//     );
//     Field<vector> velocities(particles.size(),vector::zero);
//
//     forAll(particles, i)
//     {
//         velocities[i] = velocityDict.get<vector>(particles[i]);
//     }
    
//     const vector velocity = solidInitialiseDict_.get<vector>("velocity");
    
//     const dictionary& temperatureDict
//     (
//         solidInitialiseDict_.subDict("temperatures")
//     );
//     Field<scalar> temperatures(particles.size(),0.0);
//     
//     forAll(particles, i)
//     {
//         temperatures[i] = temperatureDict.get<scalar>(particles[i]);
//     }
    
    const scalar temperature = initialiseDict_.get<scalar>("temperature");

    numberDensities /= cloud_.nRealParticles();
     
    const auto& meshCC = cloud_.mesh().cellCentres();

    forAll(mesh_.cells(), cellI)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            cellI
        );

        for (const tetIndices& cellTetIs : cellTets)
        {
            
//             const tetIndices& cellTetIs = cellTets[tetI];

            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            forAll(particles, i)
            {
                const word& particleName(particles[i]);

                label typeId(cloud_.typeIdList().find(particleName));
                
                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << particleName << "not defined." << nl
                        << abort(FatalError);
                }

                
                scalar numberDensity = numberDensities[i];

                // Calculate the number of particles required
                scalar RWF = cloud_.dsmcCloudPtr()->axiRWF(meshCC[cellI]);
                scalar particlesRequired = numberDensity*tetVolume/RWF;
                
/*                
                if(cloud_.dsmcCloudReference()->axisymmetric())
                {
                    
                    const point& cC = cloud_.mesh().cellCentres()[cellI];
                    scalar radius = cC.y();
                    
                    scalar RWF = 1.0;
                    
                    RWF = 1.0 + cloud_.dsmcCloudReference()->maxRWF()*(radius/cloud_.dsmcCloudReference()->radialExtent());
                    
                    particlesRequired /= RWF;
                }*/
                
                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);
                
                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > cloud_.rndGen().sample01<scalar>()
                )
                {
                    ++nParticlesToInsert;
                }
                
                for (label pI = 0; pI < nParticlesToInsert; ++pI)
                {
                    point p = tet.randomPoint(cloud_.rndGen());
                    
                    label newParcel = 0;
                    
                    scalar RWF = cloud_.dsmcCloudPtr()->axiRWF(meshCC[cellI]);
//                     
//                     if(cloud_.dsmcCloudReference()->axisymmetric())
//                     {                      
//                         const point& cC = cloud_.mesh().cellCentres()[cellI];
//                         scalar radius = cC.y();
//                         
//                         RWF = 1.0 + cloud_.dsmcCloudReference()->maxRWF()*(radius/cloud_.dsmcCloudReference()->radialExtent());
//                     }

                    cloud_.addNewParticle
                    (
//                         mesh_,
                        p,
                        cellI,
                        typeId,
                        newParcel,
                        temperature,
                        RWF,
                        velocities[i],
                        (vector::zero),
                        (vector::zero),
                        (vector::zero)
                    );
                }
            }
        }
    }
    
//- for interphase coupling
//     label mostAbundantType(findMax(numberDensities));
//
//     const solidParticleCoupling::constantProperties& cP = cloud_.constSolidProps
//     (
//         mostAbundantType
//     );
//
//     //- if simply use "velocity", when the "velocity" is (0 0 0), then there will always be 0 DSMC reflection candidates from solid particle surfaces. To avoid this situation, user must self-define a initialRelativeSpeed
//
//    cloud_.sigmaTcRMax().primitiveFieldRef() = cP.sigmaT()*interphaseInitialRelativeSpeed*rndGenS_.sample01<scalar>();
//
//
//     cloud_.sigmaTcRMax().correctBoundaryConditions();
    
//- for particle-particle collision
//     if(cloud_.solidCollisionDetectionModel().active()==true)
//     {
//         //- in the case of "velocity" = 0, there will always be no collision candidates in particle-particle collisions. To avoid this, user must self-define a initialRelativeSpeed.
//         cloud_.cRsigmaPPMax().primitiveFieldRef() = pi*sqr(cP.d())*interparticleInitialRelativeSpeed*rndGenS_.sample01<scalar>();
//
//
//     }
//     cloud_.cRsigmaPPMax().correctBoundaryConditions();
//
//     cloud_.updateParticleVolumeFraction();

}


} // End namespace Foam

// ************************************************************************* //
