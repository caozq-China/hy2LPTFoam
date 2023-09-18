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
    Foam::MultiphaseBinaryCollisionModel

Description
    Templated DSMC particle collision class

SourceFiles
    BinaryCollisionModel.C
    BinaryCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef InterphaseCoupling_H
#define InterphaseCoupling_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"
// #include "fvMesh.H"
#include "dsmcParcel.H"
#include "solidParticleCoupling.H"
//#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class InterphaseCouplingModel Declaration
\*---------------------------------------------------------------------------*/

class InterphaseCoupling
{
// Private data

protected:

//         const fvMesh& mesh_;

        solidParticleCouplingCloud& spc_;
        
//         Random& rndGenS_;

public:

    //- Runtime type information
    TypeName("InterphaseCoupling");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InterphaseCoupling,
        dictionary,
        (
//             const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        ),
        (spc, dict)
    );


    // Static data members

    // Constructors
    //- Construct null from owner
//     InterphaseCoupling(solidParticleCouplingCloud& owner);

    InterphaseCoupling
    (
//         const polyMesh& mesh,
        solidParticleCouplingCloud& spc,
        const dictionary& dict
    );

    // Selectors

    static autoPtr<InterphaseCoupling> New
    (
//         const polyMesh& mesh,
        solidParticleCouplingCloud& spc,
        const dictionary& dict
    );

    // Destructor
    virtual ~InterphaseCoupling() = default;

    // Member Functions
    
        //virtual void initialConfiguration() = 0;

        //- Return the collision cross section * relative velocity product
        //- monatomic gas effect on solid particles
        virtual vector monatomicForceTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const = 0;
        
        virtual scalar monatomicEnergyTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const = 0;
		
        //- polyatomic gas effect on solid particles
        virtual vector polyatomicForceTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const = 0;
        
        virtual scalar polyatomicEnergyTransferToSolidParticle
        (
            const solidParticleCoupling& pP,
            const dsmcParcel& pG,
            const scalar& cellVolume,
            const scalar& gasRWF
        ) const = 0;
        
        //- particle effect on gas molecules
        virtual void moleculePostCollisionVelocityUpdate
        (
            const solidParticleCoupling& pP,
            dsmcParcel& pG
        ) = 0;
        
        //- check if two-way coupling is not enable
        virtual bool enableMoleculeVelocityUpdate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "BinaryCollisionModelI.H"

#endif

// ************************************************************************* //
