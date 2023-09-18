/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "constants.H"
using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const Foam::List<Foam::word>& Foam::xParticleCloud::typeIdList() const
{
    return typeIdList_;
}

inline Foam::scalar Foam::xParticleCloud::nRealParticles() const
{
    return nRealParticles_;
}

inline Foam::dsmcCloud* Foam::xParticleCloud::dsmcCloudPtr() const
{
    return dsmcCloudPtr_;
}

inline const Foam::List<Foam::xParticle::constantProperties>& Foam::xParticleCloud::constProps() const
{
    return constProps_;
}

inline const Foam::xParticle::constantProperties& Foam::xParticleCloud::constProps
(
    const label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorInFunction
            << "extraterrestrial particle constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}

inline const Foam::InteractionLists<Foam::xParticle>& Foam::xParticleCloud::il() const
{
    return il_;
}

inline const Foam::CollisionModel<Foam::xParticleCloud>& Foam::xParticleCloud::collision() const
{
    return *collisionModel_;
}

inline Foam::CollisionModel<Foam::xParticleCloud>& Foam::xParticleCloud::collision()
{
    return *collisionModel_;
}

inline Foam::Random& Foam::xParticleCloud::rndGen()
{
    return rndGen_;
}

inline void Foam::xParticleCloud::clearParticles()
{
    return IDLList<xParticle>::clear();
}

// ************************************************************************* //
