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

\*---------------------------------------------------------------------------*/

#include "mppicTimeScaleModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mppicTimeScaleModels, 0);
    defineRunTimeSelectionTable(mppicTimeScaleModels, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mppicTimeScaleModels::mppicTimeScaleModels
(
    const dictionary& dict
)
:
    alphaPacked_(dict.get<scalar>("alphaPacked")),
    e_(dict.get<scalar>("e"))
{
}


// Foam::TimeScaleModel::TimeScaleModel
// (
//     const TimeScaleModel& cm
// )
// :
//     alphaPacked_(cm.alphaPacked_),
//     e_(cm.e_)
// {
// }


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mppicTimeScaleModels> Foam::mppicTimeScaleModels::New
(
    const dictionary& dict
)
{
    word modelType
    (
        dict.get<word>("TimeScaleModel")
    );

    Info<< "Selecting mppicTimeScaleModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "time scale model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<mppicTimeScaleModels>
    (
        cstrIter()(dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mppicTimeScaleModels::~mppicTimeScaleModels()
{}


// ************************************************************************* //
