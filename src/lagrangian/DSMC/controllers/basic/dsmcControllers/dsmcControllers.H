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

Class
    Foam::dsmcControllers

Description
    Stores all the information for the dsmcControllers models defined within
    the dsmcControllersDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef dsmcControllers_H
#define dsmcControllers_H

#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"

#include "dsmcStateController.H"
#include "dsmcFluxController.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dsmcControllers Declaration
\*---------------------------------------------------------------------------*/

class dsmcControllers
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        IOdictionary dsmcControllersDict_;

        //- Number of state controllers
        label nStateControllers_;

        //- Number of flux controllers
        label nFluxControllers_;

        //- State controllers
        PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List<autoPtr<dsmcStateController>> stateControllers_;

        //- Flux controllers
        PtrList<entry> fluxControllersList_;
        List<word> fCNames_;
        List<label> fCIds_;
        List<fileName> fCFixedPathNames_;
        List<autoPtr<dsmcFluxController>> fluxControllers_;


public:

    // Constructors

        //- Null Constructor
        dsmcControllers(const Time& t, const polyMesh& mesh);

        //- Constructor for dsmcFOAM
        dsmcControllers(const Time& t, const polyMesh& mesh, dsmcCloud& cloud);

    //- Destructor
    ~dsmcControllers() = default;


    // Member Functions

        //- Initial configuration
        void initialConfig();

        //- Update time information
        void updateTimeInfo();

        //- Different control stages
        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        //- Calculate controller properties
        void calculateProps();

        //- Output controller properties
        void outputResults();


        // Access

            label nStateControllers() const;

            inline List<autoPtr<dsmcStateController>>& stateControllers();

            inline const List<autoPtr<dsmcStateController>>&
                stateControllers() const;

            inline List<autoPtr<dsmcFluxController>>& fluxControllers();

            inline const List<autoPtr<dsmcFluxController>>&
                fluxControllers() const;

            inline const List<word>& stateControllersNames() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcControllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
