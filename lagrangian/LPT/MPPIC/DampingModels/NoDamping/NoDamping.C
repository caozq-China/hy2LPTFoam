/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "NoDamping.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NoDamping, 0);

addToRunTimeSelectionTable(DampingModel, NoDamping, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NoDamping::NoDamping
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
:
    DampingModel(dict, cloud)
{}


// template<class CloudType>
// Foam::DampingModels::NoDamping<CloudType>::NoDamping
// (
//     const NoDamping<CloudType>& cm
// )
// :
//     DampingModel<CloudType>(cm)
// {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NoDamping::~NoDamping()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::NoDamping::velocityCorrection
(
    solidParcel& p,
    const scalar deltaT
) const
{
    return Zero;
}


bool Foam::NoDamping::active() const
{
    return false;
}

void Foam::NoDamping::cacheFields(const bool store){}


// ************************************************************************* //
