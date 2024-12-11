/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "deletionPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(deletionPatch, 0);
addToRunTimeSelectionTable(solidPatchBoundary, deletionPatch, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deletionPatch::deletionPatch
(
    const polyMesh& mesh,
    solidParcelCloud& cloud,
    const dictionary& dict
)
:
    solidPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    allSpecies_(false),
    typeIds_()
{
    measurePropertiesAtWall_ = false;
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::deletionPatch::initialConfiguration()
{}


void Foam::deletionPatch::calculateProperties()
{}


void Foam::deletionPatch::controlParticle
(
    solidParcel& p,
    solidParcel::trackingData& td
)
{
    td.keepParticle = false;
}


void Foam::deletionPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::deletionPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    solidPatchBoundary::updateProperties(dict);

    setProperties();
}


void Foam::deletionPatch::setProperties()
{
    propsDict_.readIfPresent("allParticleTypes", allSpecies_);

    if (!allSpecies_)
    {
        typeIds_ = spc_.getTypeIDs(propsDict_);
    }
}


// ************************************************************************* //
