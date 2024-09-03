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

#include "noSlipTsirkunovWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noSlipTsirkunovWallPatch, 0);

addToRunTimeSelectionTable
(
    solidPatchBoundary,
    noSlipTsirkunovWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noSlipTsirkunovWallPatch::noSlipTsirkunovWallPatch
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

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noSlipTsirkunovWallPatch::initialConfiguration()
{}


void Foam::noSlipTsirkunovWallPatch::calculateProperties()
{}


void Foam::noSlipTsirkunovWallPatch::controlParticle
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

    scalar beta1 = 0.5*pi-acos(mag(Un)/mag(U));
    aNormal_ = 1-(1-exp(-0.1*pow(mag(U),0.61)))*sin(beta1);

    if(Un > 0)
    {
        U -= (1.0 + aNormal_)*Un*nw;
    }

    measurePropertiesAfterControl(p);
}


void Foam::noSlipTsirkunovWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::noSlipTsirkunovWallPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    solidPatchBoundary::updateProperties(dict);
    
    setProperties();
}


void Foam::noSlipTsirkunovWallPatch::setProperties()
{
    
}


// ************************************************************************* //
