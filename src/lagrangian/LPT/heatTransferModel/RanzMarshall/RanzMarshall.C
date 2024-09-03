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

#include "RanzMarshall.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RanzMarshall, 0);
    addToRunTimeSelectionTable
    (
        HeatTransferModel,
        RanzMarshall,
        dictionary
    );
};
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::RanzMarshall::RanzMarshall
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
:
    HeatTransferModel(dict, cloud, typeName)
{}


Foam::RanzMarshall::RanzMarshall(const RanzMarshall& htm)
:
    HeatTransferModel(htm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::RanzMarshall::active() const
{
    return true;
}

Foam::scalar Foam::RanzMarshall::Nu
(
    const scalar Re,
    const scalar Pr,
    const scalar Ma
) const
{
    return 2.0 + 0.6*sqrt(Re)*cbrt(Pr);
}


// ************************************************************************* //
