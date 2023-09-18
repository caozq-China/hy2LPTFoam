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

#include "solidSpecularWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidSpecularWallPatch, 0);

addToRunTimeSelectionTable(solidPatchBoundary, solidSpecularWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidSpecularWallPatch::solidSpecularWallPatch
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
:
    solidPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    calculateHeatConduction_ = false;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void solidSpecularWallPatch::initialConfiguration()
{}

void solidSpecularWallPatch::calculateProperties()
{}

void solidSpecularWallPatch::controlParticle
(
    solidParticleCoupling& p, 
    solidParticleCoupling::trackingData& td
)
{

    measurePropertiesBeforeControl(p);

    specularReflection(p);
    
    measurePropertiesAfterControl(p, 0.0);
}

void solidSpecularWallPatch::injectParticlesFromWall()
{

}

void solidSpecularWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void solidSpecularWallPatch::updateProperties(const dictionary& dict)
{
    //- the main properties should be updated first
    solidPatchBoundary::updateProperties(dict);
}


} // End namespace Foam

// ************************************************************************* //
