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
    Foam::dsmcFluxController

Description
    Basic/abstract class of a flux controller

SourceFiles
    dsmcFluxControllerI.H
    dsmcFluxController.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcFluxController_H
#define dsmcFluxController_H

#include "dsmcControllerBase.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcFluxController Declaration
\*---------------------------------------------------------------------------*/

class dsmcFluxController
:
    public dsmcControllerBase
{
protected:

    // Protected data

        scalar zoneSurfaceArea_;
        labelList internalFaces_;
        labelList processorFaces_;


    // Private Member Functions

        void setFacesInfo();


public:

    //- Runtime type information
    TypeName("dsmcFluxController");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcFluxController,
        dictionary,
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    dsmcFluxController
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    //- Selector
    static autoPtr<dsmcFluxController> New
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dsmcFluxController() = default;


    // Member Functions

        //- Return the control zone faces
        virtual const labelList& controlZone() const;

        virtual label isFaceOnControlZone(const label facei) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
