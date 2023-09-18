/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "solidSimpleDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidSimpleDiffuseWallPatch, 0);

addToRunTimeSelectionTable(solidPatchBoundary, solidSimpleDiffuseWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidSimpleDiffuseWallPatch::solidSimpleDiffuseWallPatch
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
:
    solidPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    coeffRestitution_(0.0),
    coeffFriction_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void solidSimpleDiffuseWallPatch::initialConfiguration()
{}

void solidSimpleDiffuseWallPatch::calculateProperties()
{}

void solidSimpleDiffuseWallPatch::controlParticle
(
    solidParticleCoupling& p, 
    solidParticleCoupling::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    vector nw = p.normal();
    nw.normalise();

    // Normal velocity magnitude before collision
    scalar U_dot_nw = U & nw;
    
    // Wall tangential velocity (flow direction) before collision
    vector Ut = U - U_dot_nw*nw;
    
    Ut = Ut*coeffFriction_;
    
    scalar U_dot_nw_New = U_dot_nw*coeffRestitution_;
    
    vector UnNew = -U_dot_nw_New*nw;
    
    U = Ut + UnNew;

    measurePropertiesAfterControl(p, 0.0);
}

void solidSimpleDiffuseWallPatch::injectParticlesFromWall()
{}

void solidSimpleDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void solidSimpleDiffuseWallPatch::updateProperties(const dictionary& dict)
{
    //- the main properties should be updated first
    solidPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");

    setProperties();
}

void solidSimpleDiffuseWallPatch::setProperties()
{
    coeffRestitution_ = propsDict_.get<scalar>("CoefficientOfRestitution");
    coeffFriction_ = propsDict_.get<scalar>("CoefficientOfFriction");
}
} // End namespace Foam

// ************************************************************************* //
