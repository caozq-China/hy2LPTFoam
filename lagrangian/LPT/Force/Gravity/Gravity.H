/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::Gravity

Group

Description
    Calculates particle gravity force

SourceFiles
    GravityI.H
    Gravity.C

\*---------------------------------------------------------------------------*/

#ifndef Gravity_H
#define Gravity_H

#include "Force.H"
#include "solidParcel.H"
#include "solidParcelCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class Gravity Declaration
\*---------------------------------------------------------------------------*/
class Gravity
:
    public Force
{
    // Private data

        //- Reference to the acceleration due to gravity
        const vector& g_;


public:

    //- Runtime type information
    TypeName("gravity");


    // Constructors

        //- Construct from mesh
        Gravity
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        Gravity(const Gravity& gf);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new Gravity(*this)
            );
        }


    //- Destructor
    virtual ~Gravity();


    // Member Functions

        // Access

            //- Return the acceleration due to gravity
            inline const vector& g() const;


        // Evaluation
            //- Calculate the non-coupled force
            forceSuSp calcNonCoupled
            (
                const solidParcel& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
