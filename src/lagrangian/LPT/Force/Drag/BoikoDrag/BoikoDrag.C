/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "BoikoDrag.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BoikoDrag, 0);
    addToRunTimeSelectionTable
    (
        Force,
        BoikoDrag,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::BoikoDrag::Cd
(
    const scalar Re,
    const scalar Ma12
) const
{

    return (0.38 +24/Re + 4/sqrt(Re))*(1+exp(-0.43/pow(Ma12, 4.67)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BoikoDrag::BoikoDrag
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, true)
{}


Foam::BoikoDrag::BoikoDrag
(
    const BoikoDrag& df
)
:
    Force(df)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BoikoDrag::~BoikoDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::forceSuSp Foam::BoikoDrag::calcCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const scalar speedOfSound = mag(p.Uc())/p.Mac();
    const scalar Ma12 = mag(p.U()-p.Uc())/speedOfSound;
    //- Re is based on gas-particle relative velocity
    scalar cd = Cd(Re,Ma12);
 
    value.Sp() = (mass/p.rho())*0.75*cd*Re*muc/sqr(p.d());
    return value;
}

// ************************************************************************* //
