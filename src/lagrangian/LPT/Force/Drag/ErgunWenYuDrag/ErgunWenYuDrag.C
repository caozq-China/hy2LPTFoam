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

#include "ErgunWenYuDrag.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ErgunWenYuDrag, 0);
    addToRunTimeSelectionTable
    (
        Force,
        ErgunWenYuDrag,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::ErgunWenYuDrag::CdRe
(
    const scalar Re
) const
{

    if (Re > 1000.0)
    {
        return 0.44*Re;
    }
    else
    {
        return 24.0*(1.0 + 0.15*pow(Re, 0.687));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ErgunWenYuDrag::ErgunWenYuDrag
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, true),
    thetaC_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            coeffDict_.lookup("thetaC")
        )
    )
{}


Foam::ErgunWenYuDrag::ErgunWenYuDrag
(
    const ErgunWenYuDrag& df
)
:
    Force(df),
    thetaC_
    (
        this->mesh().template lookupObject<volScalarField>
        (
            coeffDict_.lookup("thetaC")
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ErgunWenYuDrag::~ErgunWenYuDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::forceSuSp Foam::ErgunWenYuDrag::calcCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    scalar thetaC(thetaC_[p.cell()]);

    if (thetaC < 0.8)
    {
        value.Sp() = (mass/p.rho())
                    *(150.0*(1.0 - thetaC)/thetaC + 1.75*Re)*muc/(thetaC*sqr(p.d()));
    }
    else
    {
        value.Sp() = (mass/p.rho())
                    *0.75*CdRe(thetaC*Re)*muc*pow(thetaC, -2.65)/(thetaC*sqr(p.d()));
    }
 
    return value;
}

// ************************************************************************* //
