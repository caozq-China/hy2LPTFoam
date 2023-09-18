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
    dsmcControllers

Description

	Stores all the information for the dsmcControllers models defined within 
    the dsmcControllersDict, and selects & builds the models automatically.

    

\*---------------------------------------------------------------------------*/

#ifndef solidControllers_H
#define solidControllers_H


#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"

#include "solidStateController.H"
// #include "dsmcFluxController.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcControllers Declaration
\*---------------------------------------------------------------------------*/

class solidControllers
{

    // Private data

        //- The entire dictionary (containing multiple subDictionaries)


        const Time& time_;

        IOdictionary solidControllersDict_;

        // number of dsmcControllers
        label nStateControllers_;
//         label nFluxControllers_;

        //- state dsmcControllers
        PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List< autoPtr<solidStateController> > stateControllers_;

        //- flux dsmcControllers
//         PtrList<entry> fluxControllersList_;
//         List<word> fCNames_;
//         List<label> fCIds_;
//         List<fileName> fCFixedPathNames_;
//         List< autoPtr<dsmcFluxController> > fluxControllers_;

        
public:

    // Constructors

        //- Null Constructor 
        solidControllers
        (
            const Time& t,
            const polyMesh& mesh
        );

        //- Constructor for dsmcFOAM
        solidControllers
        (
            const Time& t,
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud
        );


        ~solidControllers() = default;

    // Member Functions

        //- initial configuration
        void initialConfig();

        //- update time information
        void updateTimeInfo();

        //- different control stages 
        void controlBeforeMove();
        void controlAfterMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        //- calculate controller properties
        void calculateProps();

        //- output controller properties
        void outputResults();
        
        // Access

            label nStateControllers() const;

            inline List< autoPtr<solidStateController> >& stateControllers();
            
            inline const List< autoPtr<solidStateController> >& stateControllers() const;

//             inline List< autoPtr<dsmcFluxController> >& fluxControllers();
//             inline const List< autoPtr<dsmcFluxController> >& fluxControllers() const;

            inline const List<word>& stateControllersNames() const; 
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidControllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
