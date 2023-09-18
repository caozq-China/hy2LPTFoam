/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "NoPacking.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// class solidParticleCouplingCloud;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NoPacking, 0);

addToRunTimeSelectionTable(mppicPackingModel, NoPacking, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NoPacking::NoPacking
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    mppicPackingModel(spc,dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::NoPacking::~NoPacking()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::vector Foam::NoPacking::velocityUpdateAndCorrection
(
    solidParticleCoupling& p,
    const scalar deltaT
) const
{
    return Zero;
}


bool Foam::NoPacking::active() const
{
    return false;
}

}

void Foam::NoPacking::cacheFields(const bool store)
{}
// ************************************************************************* //
