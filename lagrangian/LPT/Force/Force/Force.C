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

\*---------------------------------------------------------------------------*/

#include "Force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Force, 0);

    defineRunTimeSelectionTable(Force, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Force::Force
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceTypeName,
    const bool readCoeffs
)
:
    cloud_(cloud),
    mesh_(mesh),
    coeffDict_(dict)
{
}

Foam::Force::Force(const Force& pf)
:
    cloud_(pf.cloud_),
    mesh_(pf.mesh_),
    coeffDict_(pf.coeffDict_)
{}

Foam::autoPtr<Foam::Force>
Foam::Force::New
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict,//particleProperties_.subDict(Forces)
    const word& forceTypeName
)
{
    Info<< "    Selecting particle force " << forceTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(forceTypeName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown particle force type "
            << forceTypeName
            << ", constructor not in hash table" << nl << nl
            << "    Valid particle force types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Force>
    (
        cstrIter()
        (
            cloud,
            mesh,
            dict
        )
    );

}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
Foam::Force::~Force()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::Force::cacheFields(const bool store)
{}

Foam::forceSuSp Foam::Force::calcCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value;
    value.Su() = Zero;
    value.Sp() = 0.0;
    return value;
}

Foam::forceSuSp Foam::Force::calcNonCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value;
    value.Su() = Zero;
    value.Sp() = 0.0;
    return value;
}

Foam::scalar Foam::Force::massAdd
(
    const solidParcel& p,
    const scalar mass
) const
{
    return 0.0;
}


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //从

// ************************************************************************* //
