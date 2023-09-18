/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::VariableHardSphere

Description
    Variable Hard Sphere BinaryCollision Model

\*---------------------------------------------------------------------------*/

#ifndef stochasticHardSphere_H
#define stochasticHardSphere_H

#include "solidBinaryCollisionModel.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class solidHardSphereModel Declaration
\*---------------------------------------------------------------------------*/


class stochasticHardSphere
:
    public solidBinaryCollisionModel
{
    // Private data

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- resituation coefficient
        const scalar e_;
        
        //- friction coefficient
        //const scalar f_;
        
        //const bool enableParticleSlide_;
        

public:

    //- Runtime type information
    TypeName("stochasticHardSphere");


    // Constructors

        //- Construct from dictionary
        stochasticHardSphere
        (
            const dictionary& dict,
            solidParticleCouplingCloud& cloud
        );


    //- Destructor
    virtual ~stochasticHardSphere() = default;


    // Member Functions

        //- Flag to indicate whether model activates collision model
        bool active() const;

        //- Return the collision cross section * relative velocity product
//         scalar sigmaTcR
//         (
//             const solidParticleCoupling& pP,
//             const solidParticleCoupling& pQ
//         ) const;

        //- Apply collision
        virtual void collide
        (
            solidParticleCoupling& pP,
            solidParticleCoupling& pQ
        );

        virtual void velocityCorrection
        (
            int nCorrectionStep,
            solidParticleCoupling& pP,
            solidParticleCoupling& pQ
        );

        const dictionary& coeffDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "VariableHardSphere.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
