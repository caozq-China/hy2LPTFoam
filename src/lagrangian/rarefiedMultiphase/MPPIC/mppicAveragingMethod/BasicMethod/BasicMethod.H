/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::AveragingMethods::Basic

Description
    Basic lagrangian averaging procedure.

    This is a cell-volume based average. Point values are summed over the
    computational cells and the result is divided by the cell volume.
    
    Interpolation is done assuming a constant value over each cells. Cell
    gradients are calculated by the default fvc::grad scheme, and are also
    assumed constant when interpolated.

SourceFiles
    Basic.C

\*---------------------------------------------------------------------------*/

#ifndef BasicMethod_H
#define BasicMethod_H

#include "mppicAveragingMethod.H"
#include "pointMesh.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mppicAveragingMethods
{

/*---------------------------------------------------------------------------*\
                        Class BasicMethod Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BasicMethod
:
    public mppicAveragingMethod<Type>
{
public:

    //- Public typedefs

        //- Gradient type
        typedef typename mppicAveragingMethod<Type>::TypeGrad TypeGrad;


private:

    //- Private data

        //- Cell average field
        Field<Type>& data_;

        //- Gradient field
        Field<TypeGrad> dataGrad_;


    //- Private member functions

        //- Re-calculate gradient
        virtual void updateGrad();


public:

    //- Runtime type information
    TypeName("BasicMethod");


    //- Constructors

        //- Construct from components
        BasicMethod
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        BasicMethod(const BasicMethod<Type>& am);

        //- Construct and return a clone
        virtual autoPtr<mppicAveragingMethod<Type> > clone() const
        {
            return autoPtr<mppicAveragingMethod<Type> >
            (
                new BasicMethod<Type>(*this)
            );
        }


    //- Destructor
    virtual ~BasicMethod();


    //- Member Functions

        //- Add point value to interpolation
        void add
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const Type& value
        );

        //- Interpolate
        Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;

        //- Interpolate gradient
        TypeGrad interpolateGrad
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;

        //- Return an internal field of the average
        tmp<Field<Type> > primitiveField() const;

        //- Return an internal field of the gradient
        tmp<Field<TypeGrad> > internalFieldGrad() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mppicAveragingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BasicMethod.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
