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

#include "zoneFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(zoneFill, 0);

addToRunTimeSelectionTable(configuration, zoneFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
zoneFill::zoneFill
(
    solidParcelCloud& cloud,
    const dictionary& dict
)
:
    configuration(cloud, dict),
    randomNormalizedVelocity_(Switch(dict.lookup("randomNormalizedVelocity")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void zoneFill::setInitialConfiguration()
{

    Info<< nl << "Initialising particles" << endl;

    const vector velocity(initialiseDict_.lookup("velocity"));

    const dictionary& numberDensitiesDict
    (
        initialiseDict_.subDict("numberDensities")
    );

    List<word> particles(numberDensitiesDict.toc());

    Field<scalar> numberDensities(particles.size());

    forAll(particles, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(particles[i])
        );
    }

    numberDensities /= cloud_.nParticle();

    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word zoneName(initialiseDict_.lookup("zoneName"));
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

        forAll(zone,c)
        {
            const label& cellI = zone[c];
    
            List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
            (
                mesh_,
                cellI
            );
    
            forAll(cellTets, tetI)
            {
                const tetIndices& cellTetIs = cellTets[tetI];
    
                tetPointRef tet = cellTetIs.tet(mesh_);
    
                scalar tetVolume = tet.mag();
    
                forAll(particles, i)
                {
                    const word& particleName(particles[i]);
    
                    label typeId(findIndex(cloud_.typeIdList(),particleName));
    
                    if (typeId == -1)
                    {
                        FatalErrorInFunction
                            << "typeId " << particleName << "not defined." << nl
                            << abort(FatalError);
                    }
    
                    scalar numberDensity = numberDensities[i];
                    
                    scalar particlesRequired = numberDensity*tetVolume;

                    if(cloud_.axisymmetric())
                    {                    
                        const point& cC = cloud_.mesh().cellCentres()[cellI];
                        scalar radius = cC.y();
                        
                        scalar RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                        particlesRequired /= RWF;
                    }

                    //particlesRequired /= cloud_.nParticle(cellI);

                    label nParticlesToInsert = label(particlesRequired);

                    // Add another particle with a probability proportional to the
                    // remainder of taking the integer part of particlesRequired
                    if
                    (
                        (particlesRequired - nParticlesToInsert)
                            > cloud_.rndGen().sample01<scalar>()
                    )
                    {
                        nParticlesToInsert++;
                    }
    
                    for (label pI = 0; pI < nParticlesToInsert; pI++)
                    {
                        point position = tet.randomPoint(cloud_.rndGen());

                        vector U = velocity;
                        if(randomNormalizedVelocity_)
                        {  
                            
                            scalar epsilonAngle;
                            epsilonAngle = cloud_.rndGen().sample01<scalar>() * 2.0 * pi;
                            scalar cosElevationAngle = 2.0*cloud_.rndGen().sample01<scalar>()-1;
                            scalar sinElevationAngle = sqrt(1-sqr(cosElevationAngle));

                            U.x() = cosElevationAngle;
                            U.y() = sinElevationAngle*cos(epsilonAngle);
                            U.z() = sinElevationAngle*sin(epsilonAngle);
                        }

                        //const scalar& RWF = cloud_.coordSystem().RWF(cellI);

                        scalar RWF = 1.0;
                        
                        if(cloud_.axisymmetric())
                        {                      
                            const point& cC = cloud_.mesh().cellCentres()[cellI];
                            scalar radius = cC.y();
                            
                            RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        }

                        cloud_.addNewParcel
                        (
                            mesh_,
                            cloud_.constProps(typeId),
                            position,
                            U,
                            RWF,
                            cellI,
                            cellTetIs.face(),
                            cellTetIs.tetPt(),
                            typeId,
                            -1    
                        );
                    }
                }
            }
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
