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

#include "Fox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Fox, 0);
    addToRunTimeSelectionTable
    (
        HeatTransferModel,
        Fox,
        dictionary
    );
};
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Fox::Fox
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
:
    HeatTransferModel(dict, cloud, typeName)
{}


Foam::Fox::Fox(const Fox& htm)
:
    HeatTransferModel(htm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::Fox::~Fox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::Fox::active() const
{
    return true;
}

Foam::scalar Foam::Fox::Nu
(
    const scalar Re,
    const scalar Pr,
    const scalar Ma
) const
{
    return 2*exp(-Ma)/(1+17*Ma/Re)+0.459*pow(Pr,0.33)*pow(Re,0.55)*(1+0.5*exp(-17*Ma/Re))/1.5;
}


// ************************************************************************* //
