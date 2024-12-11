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

#include "reboundWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(reboundWallPatch, 0);

addToRunTimeSelectionTable
(
    solidPatchBoundary,
    reboundWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reboundWallPatch::reboundWallPatch
(
    const polyMesh& mesh,
    solidParcelCloud& cloud,
    const dictionary& dict
)
:
    solidPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    aNormal_(0.0),
    aTangential_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    propsDict_ = dict.subDict(typeName + "Properties");
    aNormal_ = readScalar(propsDict_.lookup("normalCoefficientOfRestitution"));
    aTangential_ = readScalar(propsDict_.lookup("tangentialCoefficientOfRestitution"));

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reboundWallPatch::initialConfiguration()
{}


void Foam::reboundWallPatch::calculateProperties()
{}


void Foam::reboundWallPatch::controlParticle
(
    solidParcel& p,
    solidParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();
    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude before collision
    scalar Un = U & nw;
    
    // Wall tangential velocity (flow direction) before collision
    vector Ut = U - Un*nw;

    if(Un > 0)
    {
        U -= (1.0 + aNormal_)*Un*nw;
    }

    U -= (1-aTangential_)*Ut;

    measurePropertiesAfterControl(p);
}


void Foam::reboundWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::reboundWallPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    solidPatchBoundary::updateProperties(dict);

    setProperties();
}


void Foam::reboundWallPatch::setProperties()
{
    
}


// ************************************************************************* //
