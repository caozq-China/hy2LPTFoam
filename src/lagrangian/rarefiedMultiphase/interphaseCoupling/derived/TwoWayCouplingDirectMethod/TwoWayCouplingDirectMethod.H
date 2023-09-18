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

#ifndef TwoWayCouplingDirectMethod_H
#define TwoWayCouplingDirectMethod_H

#include "InterphaseCoupling.H"
#include "solidParticleCouplingCloud.H"
//#include "dsmcCloud.H"
//#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class OneWayCoupling Declaration
\*---------------------------------------------------------------------------*/


class TwoWayCouplingDirectMethod
:
    public InterphaseCoupling
{
public:

    //- Runtime type information
    TypeName("TwoWayCouplingDirectMethod");


    // Constructors

        //- Construct from dictionary
        TwoWayCouplingDirectMethod
        (
//             const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );


    //- Destructor
    virtual ~TwoWayCouplingDirectMethod() = default;


    // Member Functions
    
        //void initialConfiguration()
        //- monatomic gas effect on solid particles
        vector monatomicForceTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const;
        
        scalar monatomicEnergyTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const;
        
        vector polyatomicForceTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const;
        
        scalar polyatomicEnergyTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const;
        
        void moleculePostCollisionVelocityUpdate
        (
            const solidParticleCoupling& pP,
            dsmcParcel& pG
        );
        
        bool enableMoleculeVelocityUpdate();
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
