/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HeatTransferModel, 0);

    defineRunTimeSelectionTable(HeatTransferModel, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::HeatTransferModel::HeatTransferModel(solidParcelCloud& cloud)
:
    BirdCorrection_(false)
{}


Foam::HeatTransferModel::HeatTransferModel
(
    const dictionary& dict,
    solidParcelCloud& cloud,
    const word& type
)
:
    BirdCorrection_(dict.subDict(typeName+"Properties").lookup("BirdCorrection"))
{}


Foam::HeatTransferModel::HeatTransferModel
(
    const HeatTransferModel& htm
)
:
    BirdCorrection_(htm.BirdCorrection_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HeatTransferModel::~HeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Switch& Foam::HeatTransferModel::BirdCorrection() const
{
    return BirdCorrection_;
}

Foam::scalar Foam::HeatTransferModel::htc
(
    const scalar dp,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW,
    const scalar Ma
) const
{
    const scalar Nu = this->Nu(Re, Pr, Ma);

    scalar htc = Nu*kappa/dp;

    if (BirdCorrection_ && (mag(htc) > ROOTVSMALL) && (mag(NCpW) > ROOTVSMALL))
    {
        const scalar phit = min(NCpW/htc, 50);
        if (phit > 0.001)
        {
            htc *= phit/(exp(phit) - 1.0);
        }
    }

    return htc;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HeatTransferModelNew.C"

// ************************************************************************* //
