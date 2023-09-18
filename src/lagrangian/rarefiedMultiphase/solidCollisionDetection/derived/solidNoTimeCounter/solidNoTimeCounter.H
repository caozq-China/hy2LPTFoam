/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
 

Description

    no-time counter with reactions

SourceFiles
    noTimeCounterI.H
    noTimeCounter.C
    noTimeCounterIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidNoTimeCounter_H
#define solidNoTimeCounter_H


#include "solidCollisionDetection.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidNoTimeCounter Declaration
\*---------------------------------------------------------------------------*/

class solidNoTimeCounter
:
    public solidCollisionDetection
{

private:

    // Private data

//         dictionary propsDict_;
    label infoCounter_;




public:

    //- Runtime type information
    TypeName("solidNoTimeCounter");

    // Static data members

    // Constructors

        //- construct from components
        solidNoTimeCounter
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        virtual ~solidNoTimeCounter() = default;

    // Member Functions
            bool active() const;

        // Access

        // Check

        // Edit

            void initialConfiguration();

            void collide();

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
