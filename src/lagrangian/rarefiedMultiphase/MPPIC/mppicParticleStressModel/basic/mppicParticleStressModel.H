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
    An abstract class for selecting collisions in a DSMC simulation

Description

SourceFiles
    collisionPartnerSelectionI.H
    collisionPartnerSelection.C
    collisionPartnerSelectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef mppicParticleStressModel_H
#define mppicParticleStressModel_H

#include "Field.H"
#include "FieldField.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "solidParticleCoupling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class collisionPartnerSelection Declaration
\*---------------------------------------------------------------------------*/

class mppicParticleStressModel
{
private:
    
    //- No copy assignment
        void operator=(const mppicParticleStressModel&) = delete;
    
    
protected:
    
    //- Close pack volume fraction
    scalar alphaPacked_;

public:

    //- Runtime type information
    TypeName("mppicParticleStressModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            mppicParticleStressModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Static data members

    // Constructors

        mppicParticleStressModel
        (
            const dictionary& dict
        );


    // Selectors

        static autoPtr<mppicParticleStressModel> New
        (
            const dictionary& dict
        );

    // Destructor

        virtual ~mppicParticleStressModel();

    // Member Functions
        
        //- Access max volume fraction
        scalar alphaPacked() const;
        
        //- Collision stress
        virtual tmp<Field<scalar> > tau
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const = 0;

        //- Collision stress derivaive w.r.t. the volume fraction
        virtual tmp<Field<scalar> > dTaudTheta
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const = 0;

        //- Collision stress using FieldFields
        tmp<FieldField<Field, scalar> > tau
        (
            const FieldField<Field, scalar>& alpha,
            const FieldField<Field, scalar>& rho,
            const FieldField<Field, scalar>& uRms
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "collisionPartnerSelectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
