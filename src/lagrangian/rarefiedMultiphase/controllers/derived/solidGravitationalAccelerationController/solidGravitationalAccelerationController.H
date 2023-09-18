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

\*---------------------------------------------------------------------------*/

#ifndef solidGravitationalAccelerationController_H
#define solidGravitationalAccelerationController_H

#include "solidStateController.H"
#include "vector.H"
// #include "fvMesh.H"
// #include "volFields.H"
#include "solidParticleCouplingCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class solidGravitationalAccelerationController Declaration
\*---------------------------------------------------------------------------*/

class solidGravitationalAccelerationController
:
    public solidStateController
{

        //- subDictionary containing the properties
        dictionary propsDict_;

        vector acceleration_;

        // Private Member Functions
        void setProperties();

public:

        //- Runtime type information
        TypeName("solidGravitationalAccelerationController");

        // Constructors
        //- Construct from components
        solidGravitationalAccelerationController
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud,
            const dictionary& dict
        );

        // Destructor
        virtual ~solidGravitationalAccelerationController() = default;

        // Member Functions
        virtual void initialConfiguration();
        virtual void controlSolidParticlesBeforeMove();
        virtual void controlSolidParticlesAfterMove();
        virtual void controlSolidParticlesBeforeCollisions();
        virtual void controlSolidParticlesAfterCollisions();

        //- calculate any required properties
        virtual void calculateProperties();

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&);
};

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
