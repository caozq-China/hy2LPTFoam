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
    dsmcStateController

Description
    Basic/abstract class of a state controller


SourceFiles
    dsmcStateControllerI.H
    dsmcStateController.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcStateController_H
#define dsmcStateController_H

#include "dsmcControllerBase.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcStateController Declaration
\*---------------------------------------------------------------------------*/

class dsmcStateController
:
    public dsmcControllerBase
{

protected:

    // Protected data

        scalar timePeriod_;
        scalar initialTime_;


public:

    //- Runtime type information
    TypeName("dsmcStateController");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcStateController,
        dictionary,
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    dsmcStateController
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );


    //- Selector
    static autoPtr<dsmcStateController> New
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~dsmcStateController() = default;


    // Member Functions

        virtual void controlParcelsBeforeMove() = 0;
        virtual void controlParcelsBeforeCollisions() = 0;
        virtual void controlParcelsAfterCollisions() = 0;

        virtual void updateTime();

        virtual void updateProperties(const dictionary& dict);

        //- Return the control zone cells
        virtual const labelList& controlZone() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
