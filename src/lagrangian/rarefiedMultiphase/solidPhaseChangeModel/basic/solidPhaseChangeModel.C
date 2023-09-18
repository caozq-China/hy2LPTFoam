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
    collisionPartnerSelection

Description

\*----------------------------------------------------------------------------*/

#include "solidPhaseChangeModel.H"
// #include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidPhaseChangeModel, 0);

defineRunTimeSelectionTable(solidPhaseChangeModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
solidPhaseChangeModel::solidPhaseChangeModel(solidParticleCouplingCloud& owner)
:
//     dict_(dictionary::null),
    spc_(owner)
{}

solidPhaseChangeModel::solidPhaseChangeModel
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
//     mesh_(refCast<const fvMesh>(mesh)),
    spc_(spc)
    //rndGenS_(spc_.rndGenS())
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidPhaseChangeModel> solidPhaseChangeModel::New
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
{
    const word modelType
    (
        dict.getOrDefault<word>(typeName,"NoPhaseChange")
    );

    Info<< "Selecting solidPhaseChangeModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);
    
    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidPhaseChangeModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    
    return autoPtr<solidPhaseChangeModel>
    (
        cstrIter()(spc, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
