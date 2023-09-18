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

#include "NoMppicDamping.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// class solidParticleCouplingCloud;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NoMppicDamping, 0);

addToRunTimeSelectionTable(mppicDampingModel, NoMppicDamping, dictionary);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NoMppicDamping::NoMppicDamping
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    mppicDampingModel(spc,dict)
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

Foam::NoMppicDamping::~NoMppicDamping()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::NoMppicDamping::velocityCorrection
(
    solidParticleCoupling& p,
    const scalar deltaT
) const
{
    return Zero;
}


bool Foam::NoMppicDamping::active() const
{
    return false;
}

void Foam::NoMppicDamping::cacheFields(const bool store)
{
    
}
// ************************************************************************* //
