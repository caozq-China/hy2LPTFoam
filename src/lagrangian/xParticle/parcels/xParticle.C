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

// #include "xParticle.H"
#include "xParticleCloud.H"
// #include "meshTools.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<xParticle>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::xParticle::move
(
    xParticleCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{

    td.switchProcessor = false;
    td.keepParticle = true;

    const constantProperties& constProps(cloud.constProps(typeID_));

    if(td.part()==0)
    {
        // leapfrog velocity adjust part, half time-step
        scalar mass = constProps.mass();
        U_ += 0.5*trackTime*f_/mass;
        angularMomentum_ += 0.5*trackTime*torque_;
    }
    else if (td.part() == 1)
    {
        // Leapfrog tracking part
        while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
        {
            const scalar f = 1 - stepFraction();
            trackToAndHitFace(f*trackTime*U_, f, cloud, td);
        }
    }
    else
    {
        FatalErrorInFunction
            << td.part() << " is an invalid part of the integration method."
            << abort(FatalError);
    }

    return td.keepParticle;
}

bool Foam::xParticle::hitPatch
(
    xParticleCloud&,
    trackingData&
)
{
    return false;
}

void Foam::xParticle::hitProcessorPatch
(
    xParticleCloud&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::xParticle::hitWallPatch
(
    xParticleCloud& cloud,
    trackingData& td
)
{
//     const label& patchModelId = cloud.boundaries().
//     solidPatchToModelIds()[patch()];
//
//     // apply a boundary model when a molecule collides with this poly patch
//     cloud.boundaries().solidPatchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}


void Foam::xParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);

    f_ = transform(T, f_);

    angularMomentum_ = transform(T, angularMomentum_);

    torque_ = transform(T, torque_);
}


void Foam::xParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}

// ************************************************************************* //
