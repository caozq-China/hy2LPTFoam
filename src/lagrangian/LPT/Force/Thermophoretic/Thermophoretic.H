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
    Foam::Thermophoretic

Group

Description
    Calculates particle Thermophoretic force

SourceFiles
    ThermophoreticI.H
    Thermophoretic.C

\*---------------------------------------------------------------------------*/

#ifndef Thermophoretic_H
#define Thermophoretic_H

#include "Force.H"
#include "solidParcel.H"
#include "solidParcelCloud.H"
#include "interpolation.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class Thermophoretic Declaration
\*---------------------------------------------------------------------------*/
class Thermophoretic
:
    public Force
{
    // Private data

        //- Name of Temperature field
        const word TName_;

        //- Carrier phase temperature gradient interploator
        autoPtr<interpolation<vector>> gradTInterpPtr_;


public:

    //- Runtime type information
    TypeName("Thermophoretic");


    // Constructors

        //- Construct from mesh
        Thermophoretic
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        Thermophoretic(const Thermophoretic& gf);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new Thermophoretic(*this)
            );
        }


    //- Destructor
    virtual ~Thermophoretic();


    // Member Functions

        // Access
            inline const interpolation<vector>& gradTInterp() const
            {
                if (!gradTInterpPtr_.valid())
                {
                    FatalErrorInFunction
                        << "Carrier phase gradT interpolation object not set"
                        << abort(FatalError);
                }

                return gradTInterpPtr_();
            }

            //- Cache fields
            virtual void cacheFields(const bool store);


        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
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
