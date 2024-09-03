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

#include "CliftGauvinDrag.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CliftGauvinDrag, 0);
    addToRunTimeSelectionTable
    (
        Force,
        CliftGauvinDrag,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::CliftGauvinDrag::Cd
(
    const scalar Re
) const
{

    return 24*(1+0.15*pow(Re,0.687)+0.0175*pow(Re,2.16)/(pow(Re,1.16)+42500))/Re;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CliftGauvinDrag::CliftGauvinDrag
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, true)
{}


Foam::CliftGauvinDrag::CliftGauvinDrag
(
    const CliftGauvinDrag& df
)
:
    Force(df)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CliftGauvinDrag::~CliftGauvinDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::forceSuSp Foam::CliftGauvinDrag::calcCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    //- Re is based on gas-particle relative velocity
    scalar cd = Cd(Re);
 
    value.Sp() = (mass/p.rho())*0.75*cd*Re*muc/sqr(p.d());
    return value;
}

// ************************************************************************* //
