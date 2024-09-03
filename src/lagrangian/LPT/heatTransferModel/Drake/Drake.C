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

#include "Drake.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Drake, 0);
    addToRunTimeSelectionTable
    (
        HeatTransferModel,
        Drake,
        dictionary
    );
};
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Drake::Drake
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
:
    HeatTransferModel(dict, cloud, typeName)
{}


Foam::Drake::Drake(const Drake& htm)
:
    HeatTransferModel(htm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::Drake::~Drake()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::Drake::active() const
{
    return true;
}

Foam::scalar Foam::Drake::Nu
(
    const scalar Re,
    const scalar Pr,
    const scalar Ma
) const
{
    return 2+0.459*pow(Re,0.55)*pow(Pr,0.33);
}


// ************************************************************************* //
