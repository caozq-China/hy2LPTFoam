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

#ifndef mppicAveragingMethod_H
#define mppicAveragingMethod_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "barycentric.H"
#include "runTimeSelectionTables.H"

#include "FieldField.H"//- keep it
#include "volFields.H"//- keep it
#include "pointPatchField.H"//- keep it

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class collisionPartnerSelection Declaration
\*---------------------------------------------------------------------------*/
template<class Type>
class mppicAveragingMethod
:
    public regIOobject,
    public FieldField<Field, Type>
{
protected:
    
    //- Gradient type
    typedef typename outerProduct<vector, Type>::type TypeGrad;
    
    //- Dictionary
    const dictionary& dict_;

    const fvMesh& mesh_;

    //- Update the gradient calculation
    virtual void updateGrad();


public:

    //- Runtime type information
    TypeName("mppicAveragingMethod");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            mppicAveragingMethod,
            dictionary,
            (
                const IOobject& io,
                const dictionary& dict,
                const fvMesh& mesh
//                 solidParticleCouplingCloud& spc,
                
            ),
            (io, dict, mesh)
        );

    // Static data members

    // Constructors

        mppicAveragingMethod
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh& mesh,
            const labelList& size
        );
        
        //- Construct a copy
        mppicAveragingMethod(const mppicAveragingMethod<Type>& am);
        
        //- Construct and return a clone
        virtual autoPtr<mppicAveragingMethod<Type> > clone() const = 0;

    // Selectors

        static autoPtr<mppicAveragingMethod<Type> > New
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh& mesh
        );

    // Destructor

        virtual ~mppicAveragingMethod() = default;

    // Member Functions
        
        //- Add point value to interpolation
        virtual void add
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const Type& value
        ) = 0;

        //- Interpolate
        virtual Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const = 0;

        //- Interpolate gradient
        virtual TypeGrad interpolateGrad
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const = 0;

        //- Calculate the average
        virtual void average();
        virtual void average(const mppicAveragingMethod<scalar>& weight);

        //- Dummy write
        virtual bool writeData(Ostream&) const;

        //- Write using setting from DB
        virtual bool write(const bool valid = true) const;

        //- Return an internal field of the average
        virtual tmp<Field<Type> > primitiveField() const = 0;

        //- Assign to another average
        inline void operator=(const mppicAveragingMethod<Type>& x);

        //- Assign to value
        inline void operator=(const Type& x);

        //- Assign to tmp
        inline void operator=(tmp<FieldField<Field, Type> > x);

        //- Add-equal tmp
        inline void operator+=(tmp<FieldField<Field, Type> > x);

        //- Multiply-equal tmp
        inline void operator*=(tmp<FieldField<Field, Type> > x);

        //- Divide-equal tmp
        inline void operator/=(tmp<FieldField<Field, scalar> > x);
        // Access

        // Check

        // Edit

        // Write

            //- output data

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mppicAveragingMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "mppicAveragingMethod.C"
#endif
#endif

// ************************************************************************* //
