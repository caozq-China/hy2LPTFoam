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

#include "InterphaseCoupling.H"
#include "solidParticleCouplingCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(InterphaseCoupling, 0);

defineRunTimeSelectionTable(InterphaseCoupling, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// InterphaseCoupling::InterphaseCoupling
// (
//     solidParticleCouplingCloud& owner,
// )
// :
//     spc_(owner)
// {}

InterphaseCoupling::InterphaseCoupling
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
//     mesh_(refCast<const fvMesh>(mesh)),
    spc_(spc)
//     rndGenS_(spc_.rndGenS())
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<InterphaseCoupling> InterphaseCoupling::New
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
{
    word modelType
    (
        dict.get<word>("InterphaseCouplingModel")
    );

    Info<< "Selecting InterphaseCouplingModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "InterphaseCouplingModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<InterphaseCoupling>
    (
        cstrIter()(spc, dict)
    );
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
