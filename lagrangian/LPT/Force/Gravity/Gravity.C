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

#include "Gravity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Gravity, 0);
    addToRunTimeSelectionTable
    (
        Force,
        Gravity,
        dictionary
    );
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Gravity::Gravity
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, false),
    g_(cloud_.g().value())
{
}


Foam::Gravity::Gravity(const Gravity& gf)
:
    Force(gf),
    g_(gf.g_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
Foam::Gravity::~Gravity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const Foam::vector& Foam::Gravity::g() const
{
    return g_;
}


Foam::forceSuSp Foam::Gravity::calcNonCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    value.Su() = mass*g_*(1.0 - p.rhoc()/p.rho());
    return value;
}


// ************************************************************************* //
