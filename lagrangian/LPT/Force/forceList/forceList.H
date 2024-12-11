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
    Foam::forceList

Description
    List of particle forces

SourceFiles
    forceListI.H
    forceList.C

\*---------------------------------------------------------------------------*/

#ifndef forceList_H
#define forceList_H

#include "Force.H"
#include "List.H"
#include "solidParcel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParcelCloud;
// class solidParcel;

/*---------------------------------------------------------------------------*\
                      Class forceList Declaration
\*---------------------------------------------------------------------------*/

class forceList
:
    public PtrList<Force>
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Forces dictionary
        const dictionary dict_;

        //- Calculate coupled forces flag
        bool calcCoupled_;

        //- Calculate non-coupled forces flag
        bool calcNonCoupled_;


public:

    // Constructors

        //- Construct from mesh
        forceList
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool readFields
        );


    //- Destructor
    virtual ~forceList();


    // Member Functions

        // Access

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the forces dictionary
            inline const dictionary& dict() const;

            //- Set the calcCoupled flag
            inline void setCalcCoupled(bool flag);

            //- Set the calcNonCoupled flag
            inline void setCalcNonCoupled(bool flag);

        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const solidParcel& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const solidParcel& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the effective mass
            virtual scalar massEff
            (
                const solidParcel& p,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "forceListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
