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

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
inline Foam::solidParcel::constantProperties::constantProperties()
:
    rho0_(0.0),
    d0_(0.0),
    T0_(0.0),
    TMin_(0.0),
    TMax_(VGREAT),
    Cp0_(0.0)
{}


inline Foam::solidParcel::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    rho0_(readScalar(dict.lookup("rho0"))),
    d0_(readScalar(dict.lookup("d0"))),
    T0_(readScalar(dict.lookup("T0"))),
    TMin_(dict.lookupOrDefault<scalar>("TMin", 200.0)),
    TMax_(dict.lookupOrDefault<scalar>("TMax", 5000.0)),
    K0_(readScalar(dict.lookup("K0"))),
    Cp0_(readScalar(dict.lookup("Cp0")))
{}

inline Foam::solidParcel::solidParcel
(
    const polyMesh& mesh,
    const constantProperties& constProps,
    const vector& position,
    const vector& U0,
    const scalar RWF,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const label newParcel
)
:
    particle(mesh, position, celli, tetFacei, tetPti),
    U_(U0),
    UCorrect_(Zero),
    RWF_(RWF),
    d_(constProps.d0()),
    rho_(constProps.rho0()),
    T_(constProps.T0()),
    K_(constProps.K0()),
    Cp_(constProps.Cp0()),
    typeId_(typeId),
    newParcel_(newParcel),
    active_(1),
    Uc_(Zero),
    rhoc_(0.0),
    muc_(0.0),
    Tc_(0.0),
    kappac_(0.0),
    Cpc_(0.0),
    Mac_(0.0),
    gammac_(0.0),
    omegac_(0.0)
{}

inline Foam::solidParcel::solidParcel
(
    const polyMesh& mesh,
    const constantProperties& constProps,
    const vector& position,
    const vector& U0,
    const scalar T,
    const scalar RWF,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const label newParcel
)
:
    particle(mesh, position, celli, tetFacei, tetPti),
    U_(U0),
    UCorrect_(Zero),
    RWF_(RWF),
    d_(constProps.d0()),
    rho_(constProps.rho0()),
    T_(T),
    K_(constProps.K0()),
    Cp_(constProps.Cp0()),
    typeId_(typeId),
    newParcel_(newParcel),
    active_(1),
    Uc_(Zero),
    rhoc_(0.0),
    muc_(0.0),
    Tc_(0.0),
    kappac_(0.0),
    Cpc_(0.0),
    Mac_(0.0),
    gammac_(0.0),
    omegac_(0.0)
{}

// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //
inline Foam::scalar Foam::solidParcel::constantProperties::rho0() const
{
    return rho0_;
}

inline Foam::scalar Foam::solidParcel::constantProperties::d0() const
{
    return d0_;
}

inline Foam::scalar Foam::solidParcel::constantProperties::T0() const
{
    return T0_;
}

inline Foam::scalar Foam::solidParcel::constantProperties::TMin() const
{
    return TMin_;
}

inline Foam::scalar Foam::solidParcel::constantProperties::TMax() const
{
    return TMax_;
}

inline Foam::scalar Foam::solidParcel::constantProperties::K0() const
{
    return K0_;
}

inline Foam::scalar Foam::solidParcel::constantProperties::Cp0() const
{
    return Cp0_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// - particle constructor access
inline const Foam::vector& Foam::solidParcel::U() const
{
    return U_;
}

inline const Foam::vector& Foam::solidParcel::UCorrect() const
{
    return UCorrect_;
}

inline Foam::scalar Foam::solidParcel::RWF() const
{
    return RWF_;
}

inline Foam::scalar Foam::solidParcel::d() const
{
    return d_;
}

inline Foam::scalar Foam::solidParcel::rho() const
{
    return rho_;
}

inline Foam::scalar Foam::solidParcel::T() const
{
    return T_;
}

inline Foam::scalar Foam::solidParcel::K() const
{
    return K_;
}

inline Foam::scalar Foam::solidParcel::Cp() const
{
    return Cp_;
}

inline Foam::label Foam::solidParcel::typeId() const
{
    return typeId_;
}

inline Foam::label Foam::solidParcel::active() const
{
    return active_;
}

inline Foam::label Foam::solidParcel::newParcel() const
{
    return newParcel_;
}

inline const Foam::vector& Foam::solidParcel::Uc() const
{
    return Uc_;
}

inline Foam::scalar Foam::solidParcel::rhoc() const
{
    return rhoc_;
}

inline Foam::scalar Foam::solidParcel::muc() const
{
    return muc_;
}

inline Foam::scalar Foam::solidParcel::Tc() const
{
    return Tc_;
}

inline Foam::scalar Foam::solidParcel::kappac() const
{
    return kappac_;
}

inline Foam::scalar Foam::solidParcel::Cpc() const
{
    return Cpc_;
}

inline Foam::scalar Foam::solidParcel::Mac() const
{
    return Mac_;
}

inline Foam::scalar Foam::solidParcel::gammac() const
{
    return gammac_;
}

inline Foam::scalar Foam::solidParcel::omegac() const
{
    return omegac_;
}

//- particle constructor Edit
inline Foam::vector& Foam::solidParcel::U()
{
    return U_;
}

inline Foam::vector& Foam::solidParcel::UCorrect()
{
    return UCorrect_;
}

inline Foam::scalar& Foam::solidParcel::RWF()
{
    return RWF_;
}

inline Foam::scalar& Foam::solidParcel::d()
{
    return d_;
}

inline Foam::scalar& Foam::solidParcel::rho()
{
    return rho_;
}

inline Foam::scalar& Foam::solidParcel::T()
{
    return T_;
}

inline Foam::scalar& Foam::solidParcel::K()
{
    return K_;
}

inline Foam::scalar& Foam::solidParcel::Cp()
{
    return Cp_;
}

inline Foam::label& Foam::solidParcel::typeId()
{
    return typeId_;
}

inline Foam::label& Foam::solidParcel::active()
{
    return active_;
}

inline Foam::label& Foam::solidParcel::newParcel()
{
    return newParcel_;
}

inline Foam::scalar Foam::solidParcel::massCell
(
    const label celli
) const
{
    return rhoc_*mesh().cellVolumes()[celli];
}

inline Foam::scalar Foam::solidParcel::mass() const
{
    return rho_*volume();
}

inline Foam::scalar Foam::solidParcel::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}

inline Foam::scalar Foam::solidParcel::volume() const
{
    return volume(d_);
}

inline Foam::scalar Foam::solidParcel::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}

inline Foam::scalar Foam::solidParcel::areaS() const
{
    return areaS(d_);
}

inline Foam::scalar Foam::solidParcel::areaS(const scalar d)
{
    return pi*d*d;
}

inline Foam::scalar Foam::solidParcel::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/(muc + ROOTVSMALL);
}

inline Foam::scalar Foam::solidParcel::We
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar sigma
) const
{
    return rhoc*magSqr(U - Uc_)*d/(sigma + ROOTVSMALL);
}

inline Foam::scalar Foam::solidParcel::Eo
(
    const vector& a,
    const scalar d,
    const scalar sigma
) const
{
    vector dir = U_/(mag(U_) + ROOTVSMALL);
    return mag(a & dir)*(rho_ - rhoc_)*sqr(d)/(sigma + ROOTVSMALL);
}

// ************************************************************************* //
