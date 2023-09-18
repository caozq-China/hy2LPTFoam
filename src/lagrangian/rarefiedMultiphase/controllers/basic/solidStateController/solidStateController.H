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

Class
    solidStateController

Description
    Basic/abstract class of a state controller


SourceFiles
    solidStateControllerI.H
    solidStateController.C
    solidStateControllerIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidStateController_H
#define solidStateController_H

#include "solidControllerBase.H"
// #include "IOdictionary.H"
// #include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
// #include "vector.H"
// #include "volFields.H"
// #include "Random.H"
// #include "solidParticleCoupling.H"
// #include "timeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidStateController Declaration
\*---------------------------------------------------------------------------*/

class solidStateController
:
    public solidControllerBase
{

protected:

    // Protected data
        scalar timePeriod_;
        scalar initialTime_;


public:

    //- Runtime type information
    TypeName("solidStateController");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidStateController,
            dictionary,
            (
                const polyMesh& mesh,
                solidParticleCouplingCloud& cloud,
                const dictionary& dict
            ),
            (mesh, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        solidStateController
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<solidStateController> New
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~solidStateController() = default;


    // Member Functions

        //- control dsmcParcels at different stages of the integration time-step

        virtual void controlSolidParticlesBeforeMove() = 0;
        virtual void controlSolidParticlesAfterMove() = 0;
        virtual void controlSolidParticlesBeforeCollisions() = 0;
        virtual void controlSolidParticlesAfterCollisions() = 0;
        
        virtual void updateTime();

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
        
        virtual const labelList& controlZone() const;

//         void updateStateControllerProperties(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
