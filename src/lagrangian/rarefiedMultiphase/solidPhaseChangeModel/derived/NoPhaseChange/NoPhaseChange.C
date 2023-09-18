/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "NoPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NoPhaseChange, 0);

addToRunTimeSelectionTable
(solidPhaseChangeModel, NoPhaseChange, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
NoPhaseChange::NoPhaseChange
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidPhaseChangeModel(spc, dict)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool NoPhaseChange::active() const
{
    return false;
}


void NoPhaseChange::initialConfiguration()
{

}

void NoPhaseChange::temperatureCorrection
(
    solidParticleCoupling& pSolid,
    const scalar& QdeltaC,
    const scalar& deltaT
)
{
    label typeIdSolid = pSolid.typeID();
    
    pSolid.T() += (QdeltaC * deltaT / (spc_.constSolidProps(typeIdSolid).Cp() * spc_.constSolidProps(typeIdSolid).massSphere()));
}

void NoPhaseChange::phaseCheck(solidParticleCoupling& pSolid)
{
    
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
