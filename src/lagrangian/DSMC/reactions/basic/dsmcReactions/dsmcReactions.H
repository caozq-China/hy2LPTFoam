/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    dsmcReactions

Description

    Stores all the information for the dsmcReactions models defined within
    the dsmcReactionsDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef dsmcReactions_H
#define dsmcReactions_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "dsmcReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class dsmcReactions Declaration
\*---------------------------------------------------------------------------*/

class dsmcReactions
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary chemReactDict_;

        //- Number of dsmcReactions
        label nReactions_;

        PtrList<entry> reactionsList_;
        List<word> reactionNames_;
        List<label> reactionIds_;
        List<autoPtr<dsmcReaction>> reactions_;

        //- Given two pair typeIds, one is returned the modelId number
        //  in the reactions_ data member
        List<List<label>> pairAddressing_;

        label counter_;


public:

    // Constructors

        //- Null Constructor
        dsmcReactions
        (
            const Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFOAM
        dsmcReactions
        (
            const Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


        virtual ~dsmcReactions() = default;


    // Member Functions

        //- Initial configuration
        void initialConfiguration();

        void outputData();


        // Access

            inline label returnModelId(dsmcParcel& p, dsmcParcel& q);

            inline List<autoPtr<dsmcReaction>>& reactions();
            inline const List<autoPtr<dsmcReaction>>& reactions() const;

            inline const List<List<label>>& pairModelAddressing() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcReactionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
