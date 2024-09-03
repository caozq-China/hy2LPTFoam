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

#include "StasenkoWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(StasenkoWallPatch, 0);

addToRunTimeSelectionTable
(
    solidPatchBoundary,
    StasenkoWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::StasenkoWallPatch::StasenkoWallPatch
(
    const polyMesh& mesh,
    solidParcelCloud& cloud,
    const dictionary& dict
)
:
    solidPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    aNormal_(0.0),
    aTangential_(0.0),
    Ew_(0.0),
    Gw_(0.0),
    rhoW_(0.0),
    wallYieldPoint_(0.0),
    Ep_(0.0),
    Gp_(0.0),
    parcelYieldPoint_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    propsDict_ = dict.subDict(typeName + "Properties");

    Ew_ = readScalar(propsDict_.lookup("wallYoungModulus"));
    Gw_ = readScalar(propsDict_.lookup("wallShearModulus"));
    rhoW_ = readScalar(propsDict_.lookup("wallMaterialDensity"));
    wallYieldPoint_ = readScalar(propsDict_.lookup("wallYieldPoint"));

    Ep_ = readScalar(propsDict_.lookup("particleYoungModulus"));
    Gp_ = readScalar(propsDict_.lookup("particleShearModulus"));
    parcelYieldPoint_ = readScalar(propsDict_.lookup("particleYieldPoint"));

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::StasenkoWallPatch::initialConfiguration()
{}


void Foam::StasenkoWallPatch::calculateProperties()
{}


void Foam::StasenkoWallPatch::controlParticle
(
    solidParcel& p,
    solidParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    vector nw = p.normal();

    nw /= mag(nw);// nw.normalise();

    // Normal velocity magnitude before collision
    scalar Un = U & nw;
    
    // Wall tangential velocity (flow direction) before collision
    vector Ut = U - Un*nw;
    
    scalar ccw = sqrt(Ew_*Gw_)/rhoW_;
    scalar ccp = sqrt(Ep_*Gp_)/p.rho();
    scalar c= sqrt(ccw/ccp);
    scalar UnStar = sqrt(3*((parcelYieldPoint_*wallYieldPoint_)/(wallYieldPoint_+parcelYieldPoint_))/p.rho());
    scalar aNmin = sqr(1-sqrt(mag(Un)/UnStar));
    scalar aTauMin = c/(c+0.5);
    scalar beta1 = 0.5*pi-acos(mag(Un)/mag(U));
    aNormal_ = 1-(1-aNmin)*sin(beta1);
    aTangential_ = 1-(1-aTauMin)*pow(beta1/(pi/6),1/3)*pow((0.5*pi-beta1)/(pi/3),2/3);

    if(Un > 0)
    {
        U -= (1.0 + aNormal_)*Un*nw;
    }
    
    U -= (1-aTangential_)*Ut;

    measurePropertiesAfterControl(p);
}


void Foam::StasenkoWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::StasenkoWallPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    solidPatchBoundary::updateProperties(dict);
    
    setProperties();
}


void Foam::StasenkoWallPatch::setProperties()
{
}


// ************************************************************************* //
