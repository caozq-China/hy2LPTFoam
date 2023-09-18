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

#include "solidMeshFill.H"
#include "addToRunTimeSelectionTable.H"
// #include "IFstream.H"
// #include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidMeshFill, 0);

addToRunTimeSelectionTable(solidConfiguration, solidMeshFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidMeshFill::solidMeshFill
(
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    solidConfiguration(cloud, dict)
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidMeshFill::setInitialConfiguration()
{
   
    
    Info<< nl << "Initialising solid particles" << endl;

    const scalar interphaseInitialRelativeSpeed = solidInitialiseDict_.get<scalar>("interphaseInitialRelativeSpeed");

    const scalar interparticleInitialRelativeSpeed = solidInitialiseDict_.get<scalar>("interparticleInitialRelativeSpeed");
       

    const dictionary& numberDensitiesDict
    (
        solidInitialiseDict_.subDict("numberDensities")
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
        solidInitialiseDict_.subDict("velocity")
    );
    Field<vector> velocities(particles.size(),vector::zero);
    
    forAll(particles, i)
    {
        velocities[i] = velocityDict.get<vector>(particles[i]);
    }
    
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
    
    const scalar temperature = solidInitialiseDict_.get<scalar>("temperature");
    
    //- ratio of the crystallization front radius to particle radius
    const dictionary& CzRatioDict
    (
        solidInitialiseDict_.subDict("CzRatios")
    );
    Field<scalar> CzRatios(particles.size(),0.0);
    
    if(cloud_.solidPhaseChange().active())
    {
        forAll(particles, i)
        {
            CzRatios[i] = CzRatioDict.get<scalar>(particles[i]);
        }
    }
    
    //- particle phase check
    Field<scalar> phaseStates(particles.size(),0);

    const dictionary& phaseStatesDict
    (
        solidInitialiseDict_.subDict("phaseStates")
    );
    forAll(particles, i)
    {
        phaseStates[i] = phaseStatesDict.get<label>(particles[i]);
    }
    
    //- check
    if(cloud_.solidPhaseChange().active())
    {
        forAll(particles, i)
        {
            if(temperature< cloud_.constSolidProps(i).Tf())
            {
                if(phaseStates[i] ==2)
                {
                    FatalErrorIn("solidMeshFill::setInitialConfiguration()")
                    << "    Warning! Particle phase is defined inappropriately ! "
                    << "    Please check!"
                    << nl
                    << exit(FatalError); 
                }
            }
            else if(temperature<cloud_.constSolidProps(i).Tm())
            {
                if(phaseStates[i] == 2)
                {
                    FatalErrorIn("solidMeshFill::setInitialConfiguration()")
                    << "    Warning! Particle phase is defined inappropriately ! "
                    << "    Please check!"
                    << nl
                    << exit(FatalError); 
                }
            }
            else if(temperature >= cloud_.constSolidProps(i).Tm()) 
            {
                if(phaseStates[i] != 2 || phaseStates[i] != 3)
                {
                    FatalErrorIn("solidMeshFill::setInitialConfiguration()")
                    << "    Warning! Particle phase is defined inappropriately ! "
                    << "    Please check!"
                    << nl
                    << exit(FatalError); 
                }
                
            }
        }
    }
    else
    {
        forAll(particles, i)
        {
            if(phaseStates[i]==1 || phaseStates[i]==2 || phaseStates[i]==3)
            {
                FatalErrorIn("solidMeshFill::setInitialConfiguration()")
                << "    Warning! Particle phase change model is unenabled ! "
                << "    But the particle is not in solid phase now. "
                << "    If you want to use phase change model, please enable it in spcProperties. "
                << "    The options of phase states is shown below: "
                << "    0: pure solid phase "
                << "    1: unsteady phase and the core is liquid "
                << "    2: unsteady phase and the core is solid "
                << "    3: pure liquid phase "
                << "    Otherwise, the phaseState should be 0 for the case without phase change model. "
                << nl
                << exit(FatalError); 
            }
        }
    }
    
    //- calculating instantaneous particle diameter
    Field<scalar> Dps(particles.size());
    
    if(cloud_.solidPhaseChange().active())
    {
        forAll(particles, i)
        {
            if(temperature >= cloud_.constSolidProps(i).Tm()&& CzRatios[i] > 0.0&& CzRatios[i] < 1.0)
            {   
                
                scalar rhoLiquid = cloud_.particleLiquidDensityCorrection
                                        (
                                            cloud_.materialList()[i],
                                            temperature
                                        );
                Dps[i] = pow(
                            (6.0*cloud_.constSolidProps(i).massSphere()/pi)
                            /((1.0-pow(CzRatios[i],3.0))*cloud_.constSolidProps(i).rho()
                            +pow(CzRatios[i],3.0)*rhoLiquid)
                            ,1.0/3.0
                         );
            }
            else if(phaseStates[i] ==3 && temperature >= cloud_.constSolidProps(i).Tm())
            {
                scalar rhoLiquid = cloud_.particleLiquidDensityCorrection
                                    (
                                        cloud_.materialList()[i],
                                        temperature
                                    );
            
                Dps[i] = pow((6.0*cloud_.constSolidProps(i).massSphere()/(pi*rhoLiquid)),1.0/3.0);
            }
            else
            {
                Dps[i] = cloud_.constSolidProps(i).d();
            }

            //- parcel diameter
            Dps[i] = Dps[i]*pow(cloud_.nSolidParticles(),1/3);
        }
    }
    else
    {
        forAll(particles, i)
        {
            Dps[i] = cloud_.constSolidProps(i).d();
            //- parcel diameter
            Dps[i] = Dps[i]*pow(cloud_.nSolidParticles(),1/3);
        }
    }


    
    numberDensities /= cloud_.nSolidParticles();
     
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
                const word& moleculeName(particles[i]);
                
                label typeIdSolid(cloud_.typeIdSolidList().find(moleculeName));
                
                if (typeIdSolid == -1)
                {
                    FatalErrorInFunction
                        << "typeIdSolid " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                
                scalar numberDensity = numberDensities[i];

                // Calculate the number of particles required
                scalar RWF = cloud_.dsmcCloudReference()->axiRWF(meshCC[cellI]);
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
                  > rndGenS_.sample01<scalar>()
                )
                {
                    ++nParticlesToInsert;
                }
                
                for (label pI = 0; pI < nParticlesToInsert; ++pI)
                {
                    point p = tet.randomPoint(rndGenS_);
                    
                    label newParcel = 0;
                    
                    scalar RWF = cloud_.dsmcCloudReference()->axiRWF(meshCC[cellI]);
//                     
//                     if(cloud_.dsmcCloudReference()->axisymmetric())
//                     {                      
//                         const point& cC = cloud_.mesh().cellCentres()[cellI];
//                         scalar radius = cC.y();
//                         
//                         RWF = 1.0 + cloud_.dsmcCloudReference()->maxRWF()*(radius/cloud_.dsmcCloudReference()->radialExtent());
//                     }

                    cloud_.addNewSolidParticle
                    (
//                         mesh_,
                        p,
                        velocities[i],
                        (vector::zero),
                        (vector::zero),
                        (vector::zero),
                        Dps[i],
                        temperature,
                        RWF,
                        CzRatios[i],
                        cellI,
                        -1,
                        typeIdSolid,
                        phaseStates[i],
                        newParcel
                    );
                }
            }
        }
    }
    
//- for interphase coupling
    label mostAbundantType(findMax(numberDensities));

    const solidParticleCoupling::constantProperties& cP = cloud_.constSolidProps
    (
        mostAbundantType
    );
    
    //- if simply use "velocity", when the "velocity" is (0 0 0), then there will always be 0 DSMC reflection candidates from solid particle surfaces. To avoid this situation, user must self-define a initialRelativeSpeed
    
   cloud_.sigmaTcRMax().primitiveFieldRef() = cP.sigmaT()*interphaseInitialRelativeSpeed*rndGenS_.sample01<scalar>();


    cloud_.sigmaTcRMax().correctBoundaryConditions();
    
//- for particle-particle collision
    if(cloud_.solidCollisionDetectionModel().active()==true)
    {
        //- in the case of "velocity" = 0, there will always be no collision candidates in particle-particle collisions. To avoid this, user must self-define a initialRelativeSpeed.
        cloud_.cRsigmaPPMax().primitiveFieldRef() = pi*sqr(cP.d())*interparticleInitialRelativeSpeed*rndGenS_.sample01<scalar>();

        
    }
    cloud_.cRsigmaPPMax().correctBoundaryConditions();
    
    cloud_.updateParticleVolumeFraction();

}


} // End namespace Foam

// ************************************************************************* //
