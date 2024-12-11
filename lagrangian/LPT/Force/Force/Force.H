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
    Foam::Force

Group

Description
    Abstract base class for particle forces

SourceFiles
    ForceI.H
    Force.C
    ForceNew.C

\*---------------------------------------------------------------------------*/

#ifndef Force_H
#define Force_H

#include "dictionary.H"
#include "forceSuSp.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "solidParcel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParcelCloud;

/*---------------------------------------------------------------------------*\
                        Class Force Declaration
\*---------------------------------------------------------------------------*/
class Force
{

protected:

        //- Reference to the owner cloud
        solidParcelCloud& cloud_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Force coefficients dictionary
        const dictionary& coeffDict_;


public:

    //- Runtime type information
    TypeName("Force");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Force,
        dictionary,
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (cloud, mesh, dict)
    );

    // Constructors

        //- Construct from mesh
        Force
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceTypeName,
            const bool readCoeffs
        );

        //- Construct copy
        Force(const Force& pf);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new Force(*this)
            );
        }


    //- Destructor
    virtual ~Force();


    //- Selector
    static autoPtr<Force> New
    (
        solidParcelCloud& cloud,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& forceTypeName
    );


    // Member Functions

        // Access

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the force coefficients dictionary
            inline const dictionary& coeffDict() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the force
            virtual forceSuSp calcCoupled
            (
                const solidParcel& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            virtual forceSuSp calcNonCoupled
            (
                const solidParcel& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the added mass
            virtual scalar massAdd
            (
                const solidParcel& p,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ForceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
