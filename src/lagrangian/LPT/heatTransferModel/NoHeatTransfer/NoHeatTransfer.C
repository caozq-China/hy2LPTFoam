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

#include "NoHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NoHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        HeatTransferModel,
        NoHeatTransfer,
        dictionary
    );
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::NoHeatTransfer::NoHeatTransfer
(
    const dictionary&,
    solidParcelCloud& cloud
)
:
    HeatTransferModel(cloud)
{}

Foam::NoHeatTransfer::NoHeatTransfer
(
    const NoHeatTransfer& htm
)
:
    HeatTransferModel(htm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::NoHeatTransfer::~NoHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::NoHeatTransfer::active() const
{
    return false;
}

Foam::scalar Foam::NoHeatTransfer::Nu
(
    const scalar,
    const scalar,
    const scalar
) const
{
    return 0.0;
}

Foam::scalar Foam::NoHeatTransfer::Pr() const
{
    return 1.0;
}


// ************************************************************************* //
