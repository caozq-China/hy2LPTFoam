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

\*---------------------------------------------------------------------------*/

#include "solidBinaryCollisionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(solidBinaryCollisionModel, 0);
    
    defineRunTimeSelectionTable(solidBinaryCollisionModel, dictionary);
};



Foam::solidBinaryCollisionModel::solidBinaryCollisionModel(solidParticleCouplingCloud& owner)
:
    dict_(dictionary::null),
    spc_(owner)
{}


// 
Foam::solidBinaryCollisionModel::solidBinaryCollisionModel
(
    const dictionary& dict,
    solidParticleCouplingCloud& owner
//     const word& type
)
:
    dict_(dict),
    spc_(owner)
{
    
}

Foam::autoPtr<Foam::solidBinaryCollisionModel> Foam::solidBinaryCollisionModel::New
(
    const dictionary& dict,
    solidParticleCouplingCloud& owner
)
{
    const word modelType(dict.getOrDefault<word>(typeName, "NoBinaryCollision"));

    Info<< "Selecting solidBinaryCollisionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidBinaryCollisionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<solidBinaryCollisionModel>
    (
        cstrIter()(dict, owner)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::solidBinaryCollisionModel::dict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "BinaryCollisionModelNew.C"

// ************************************************************************* //
