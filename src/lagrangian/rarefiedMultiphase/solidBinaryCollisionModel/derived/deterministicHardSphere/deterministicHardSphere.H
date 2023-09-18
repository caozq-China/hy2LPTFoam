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

#ifndef deterministicHardSphere_H
#define deterministicHardSphere_H

#include "solidBinaryCollisionModel.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class solidHardSphereModel Declaration
\*---------------------------------------------------------------------------*/


class deterministicHardSphere
:
    public solidBinaryCollisionModel
{
    // Private data

        //- The coefficients dictionary
        const dictionary coeffDict_;

//         const label nCorrectionSteps_;

        //- resituation coefficient in normal direction
        scalar eN_;

        //- resituation coefficient in tengential direction
        scalar eT_;
        
        //- friction coefficient
        const scalar f_;

        //- normal spring constant
        const scalar kn_;

        //- damping coefficient
        const scalar etaN_;
        
        //const bool enableParticleSlide_;
        

public:

    //- Runtime type information
    TypeName("deterministicHardSphere");


    // Constructors

        //- Construct from dictionary
        deterministicHardSphere
        (
            const dictionary& dict,
            solidParticleCouplingCloud& cloud
        );


    //- Destructor
    virtual ~deterministicHardSphere() = default;


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
            solidParticleCoupling& pI,
            solidParticleCoupling& pJ
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
