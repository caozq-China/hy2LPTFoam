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

#include "solidZoneFill.H"
#include "addToRunTimeSelectionTable.H"
// #include "IFstream.H"
// #include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidZoneFill, 0);

addToRunTimeSelectionTable(solidConfiguration, solidZoneFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidZoneFill::solidZoneFill
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

void solidZoneFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const vector velocity = solidInitialiseDict_.get<vector>("velocities");

    const scalar interphaseInitialRelativeSpeed = solidInitialiseDict_.get<scalar>("interphaseInitialRelativeSpeed");

//     scalar interparticleInitialRelativeSpeed = readScalar(solidInitialiseDict_.lookup("interparticleInitialRelativeSpeed"));

//     const vector gasVelocity(solidInitialiseDict_.lookup("gasVelocity"));

    const scalar temperature = solidInitialiseDict_.get<scalar>("Temperature");

    const dictionary& numberDensitiesDict
    (
        solidInitialiseDict_.subDict("numberDensities")
    );
    List<word> particles(numberDensitiesDict.toc());
//     List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(particles.size());

    forAll(particles, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(particles[i]);
    }

    numberDensities /= cloud_.nSolidParticles();

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

    const dictionary& phaseStatesDict
    (
        solidInitialiseDict_.subDict("phaseStates")
    );

    //- particle phase check
    Field<scalar> phaseStates(particles.size(),0);
    forAll(particles, i)
    {
        phaseStates[i] = phaseStatesDict.get<label>(particles[i]);
    }

    //- phase check
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
                if(phaseStates[i] ==2)
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


    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word zoneName(solidInitialiseDict_.get<word>("solidZoneName"));
    const label zoneId = cellZones.findZoneID(zoneName);

    if(zoneId == -1)
    {
        FatalErrorInFunction
            << "Cannot find a solid zone: " << zoneName << nl << "in: "
            << mesh_.time().system()/"solidInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    if (zone.size())
    {
        Info << "Lattice in zone: " << zoneName << endl;

        const auto& meshCC = cloud_.mesh().cellCentres();

        for (const label cellI : zone)
        {
//             const label& cellI = zone[c];

            List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
            (
                mesh_,
                cellI
            );

            for (const tetIndices& cellTetIs : cellTets)
            {
//                 const tetIndices& cellTetIs = cellTets[tetI];

                tetPointRef tet = cellTetIs.tet(mesh_);

                scalar tetVolume = tet.mag();

                forAll(particles, i)
                {
                    const word& particleName(particles[i]);

                    label typeIdSolid(cloud_.typeIdSolidList().find(particleName));

                    if (typeIdSolid == -1)
                    {
                        FatalErrorInFunction
                            << "typeIdSolid " << particleName << "not defined." << nl
                            << abort(FatalError);
                    }

                    scalar numberDensity = numberDensities[i];

                    scalar RWF = cloud_.dsmcCloudReference()->axiRWF(meshCC[cellI]);

                    scalar particlesRequired = numberDensity*tetVolume/RWF;

//                     Info<<"particlesRequired = "<<particlesRequired<<endl;
                    // Only integer numbers of particles can be inserted
                    label nParticlesToInsert = label(particlesRequired);

                    // Add another particle with a probability proportional to the
                    // remainder of taking the integer part of particlesRequired
                    if
                    (
                        (particlesRequired - nParticlesToInsert)
                            > rndGenS_.sample01<scalar>()
                        /*(particlesRequired - nParticlesToInsert)
                            > 0*/
                    )
                    {
                        ++nParticlesToInsert;
                    }

//                     Info<<"nParticlesToInsert = "<<nParticlesToInsert<<endl;

                    for (label pI = 0; pI < nParticlesToInsert; ++pI)
                    {
                        point p = tet.randomPoint(rndGenS_);

                        label newParcel = 0;

                        scalar RWF = cloud_.dsmcCloudReference()->axiRWF(meshCC[cellI]);

                        cloud_.addNewSolidParticle
                        (
//                             mesh_,
                            p,
                            velocity,
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

//         Info<<"forAll(zone,c) Done"<<endl;
    }

     label mostAbundantType(findMax(numberDensities));

     const solidParticleCoupling::constantProperties& cP = cloud_.constSolidProps
     (
         mostAbundantType
     );

    if(cloud_.solidCollisionDetectionModel().active()==true)
    {
        scalar interparticleInitialRelativeSpeed = solidInitialiseDict_.get<scalar>("interparticleInitialRelativeSpeed");

        for (const label cellI : zone)
        {
//             const label& cellI = zone[c];

            //- if simply use "velocity", when the "velocity" is (0 0 0), then there will always be 0 DSMC reflection candidates from solid particle surfaces. To avoid this situation, user must self-define a initialRelativeSpeed
            cloud_.sigmaTcRMax().primitiveFieldRef()[cellI] = cP.sigmaT()*interphaseInitialRelativeSpeed;

            //- for particle-particle collision
            //- in the case of "velocity" = 0, there will always be no collision candidates in particle-particle collisions. To avoid this, user must self-define a initialRelativeSpeed.
            cloud_.cRsigmaPPMax().primitiveFieldRef()[cellI] = pi*sqr(cP.d())*interparticleInitialRelativeSpeed;

        }

        cloud_.cRsigmaPPMax().correctBoundaryConditions();
    }
    else
    {
        for (const label cellI : zone)
        {
//             const label& cellI = zone[c];

            //- if simply use "velocity", when the "velocity" is (0 0 0), then there will always be 0 DSMC reflection candidates from solid particle surfaces. To avoid this situation, user must self-define a initialRelativeSpeed
            cloud_.sigmaTcRMax().primitiveFieldRef()[cellI] = cP.sigmaT()*interphaseInitialRelativeSpeed;

        }
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    cloud_.updateParticleVolumeFraction();


}


} // End namespace Foam

// ************************************************************************* //
