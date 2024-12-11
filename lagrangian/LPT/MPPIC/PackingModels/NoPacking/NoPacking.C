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

#include "NoPacking.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(NoPacking, 0);

    addToRunTimeSelectionTable(PackingModel, NoPacking, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class CloudType>
Foam::NoPacking::NoPacking
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
:
    PackingModel(dict, cloud)
{}


// template<class CloudType>
// Foam::PackingModels::NoPacking<CloudType>::NoPacking
// (
//     const NoPacking<CloudType>& cm
// )
// :
//     PackingModel<CloudType>(cm)
// {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// template<class CloudType>
Foam::NoPacking::~NoPacking()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class CloudType>
Foam::vector Foam::NoPacking::velocityCorrection
(
    // solidParcelCloud::solidParcel& p,
    solidParcel& p,
    const scalar deltaT
) const
{
    return Zero;
}


// template<class CloudType>
bool Foam::NoPacking::active() const
{
    return false;
}

void Foam::NoPacking::cacheFields(const bool store)
{}
// ************************************************************************* //
