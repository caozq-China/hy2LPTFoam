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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
inline Foam::xParticle::constantProperties::constantProperties()
:
    d_(0),
    rho_(0)
{}


inline Foam::xParticle::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    d_(dict.get<scalar>("diameter")),
    rho_(dict.get<scalar>("rho")),
    Cp_(dict.get<scalar>("specificHeat"))
{
    mass_ = constant::mathematical::pi*d_*d_*d_*rho_/6.0;
}
	


inline Foam::xParticle::xParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeID,
    const label newParticle,
    const scalar T,
    const scalar RWF,
    const vector& U,
    // const vector& omega,
    const vector& f,
    const vector& angularMomentum,
    const vector& torque
)
:
    particle(mesh, coordinates, cellI, tetFaceI, tetPtI),
    typeID_(typeID),
    newParticle_(newParticle),
    T_(T),
    RWF_(RWF),
    U_(U),
    // omega_(omega),
    f_(f),
    angularMomentum_(angularMomentum),
    torque_(torque),
    collisionRecords_()
{}

inline Foam::xParticle::xParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label typeID,
    const label newParticle,
    const scalar T,
    const scalar RWF,
    const vector& U,
    // const vector& omega,
    const vector& f,
    const vector& angularMomentum,
    const vector& torque
)
:
    particle(mesh, position, cellI),
    typeID_(typeID),
    newParticle_(newParticle),
    T_(T),
    RWF_(RWF),
    U_(U),
    // omega_(omega),
    f_(f),
    angularMomentum_(angularMomentum),
    torque_(torque),
    collisionRecords_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//- constantProperties
inline Foam::scalar Foam::xParticle::constantProperties::d() const
{
    return d_;
}

inline Foam::scalar Foam::xParticle::constantProperties::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::xParticle::constantProperties::Cp() const
{
    return Cp_;
}
/*
inline Foam::scalar Foam::xParticle::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*(0.5*effD_)*(0.5*effD_);
}*/

inline Foam::scalar Foam::xParticle::constantProperties::mass() const
{
    return mass_;
}


// - particle constructor access
inline Foam::label Foam::xParticle::typeID() const
{
    return typeID_;
}

inline Foam::label Foam::xParticle::newParticle() const
{
    return newParticle_;
}

inline Foam::scalar Foam::xParticle::T() const
{
    return T_;
}

inline Foam::scalar Foam::xParticle::RWF() const
{
    return RWF_;
}

inline const Foam::vector& Foam::xParticle::U() const
{
    return U_;
}

inline Foam::vector Foam::xParticle::omega() const
{
    const constantProperties constprop;
    return angularMomentum_/(0.1*constprop.d()*constprop.d()*constprop.mass());
}

inline const Foam::vector& Foam::xParticle::f() const
{
    return f_;
}

inline const Foam::vector& Foam::xParticle::angularMomentum() const
{
    return angularMomentum_;
}

inline const Foam::vector& Foam::xParticle::torque() const
{
    return torque_;
}

inline const Foam::collisionRecordList& Foam::xParticle::collisionRecords() const
{
    return collisionRecords_;
}

//- particle constructor Edit
inline Foam::label& Foam::xParticle::typeID()
{
    return typeID_;
}

inline Foam::label& Foam::xParticle::newParticle()
{
    return newParticle_;
}

inline Foam::scalar& Foam::xParticle::T()
{
    return T_;
}

inline Foam::scalar& Foam::xParticle::RWF()
{
    return RWF_;
}

inline Foam::vector& Foam::xParticle::U()
{
    return U_;
}

// inline Foam::vector& Foam::xParticle::omega()
// {
//     return omega_;
// }

inline Foam::vector& Foam::xParticle::f()
{
    return f_;
}

inline Foam::vector& Foam::xParticle::angularMomentum()
{
    return angularMomentum_;
}

inline Foam::vector& Foam::xParticle::torque()
{
    return torque_;
}

inline Foam::collisionRecordList& Foam::xParticle::collisionRecords()
{
    return collisionRecords_;
}
// ************************************************************************* //
