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
    solidBoundaryMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "solidBoundaryMeasurements.H"
// #include "processorPolyPatch.H"
// #include "cyclicPolyPatch.H"
// #include "wallPolyPatch.H"
#include "solidParticleCouplingCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from mesh and cloud 
solidBoundaryMeasurements::solidBoundaryMeasurements
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


//- Construct from mesh, cloud and boolean (dsmcFoam)
solidBoundaryMeasurements::solidBoundaryMeasurements
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIdSolids_(identity(cloud_.typeIdSolidList().size())),
    momentumBF_(),
    UMeanBF_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    volumeFractionBF_(),
    numberFluxBF_(),
    massFluxBF_(),
    conductiveHeatFluxBF_()
{
    
    momentumBF_.setSize(typeIdSolids_.size());
    UMeanBF_.setSize(typeIdSolids_.size());
    rhoNBF_.setSize(typeIdSolids_.size());
    rhoMBF_.setSize(typeIdSolids_.size());
    linearKEBF_.setSize(typeIdSolids_.size());
    volumeFractionBF_.setSize(typeIdSolids_.size());
    numberFluxBF_.setSize(typeIdSolids_.size());
    massFluxBF_.setSize(typeIdSolids_.size());
    conductiveHeatFluxBF_.setSize(typeIdSolids_.size());
    
    forAll(rhoNBF_, i)
    {        
        momentumBF_[i].setSize(mesh_.boundaryMesh().size());
        UMeanBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoMBF_[i].setSize(mesh_.boundaryMesh().size());
        linearKEBF_[i].setSize(mesh_.boundaryMesh().size());
        volumeFractionBF_[i].setSize(mesh_.boundaryMesh().size());
        numberFluxBF_[i].setSize(mesh_.boundaryMesh().size());
        massFluxBF_[i].setSize(mesh_.boundaryMesh().size());
        conductiveHeatFluxBF_[i].setSize(mesh_.boundaryMesh().size());
        
        forAll(rhoNBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            momentumBF_[i][j].setSize(patch.size(),vector::zero);
            UMeanBF_[i][j].setSize(patch.size(),vector::zero);
            rhoNBF_[i][j].setSize(patch.size(),0.0);
            rhoMBF_[i][j].setSize(patch.size(),0.0);
            linearKEBF_[i][j].setSize(patch.size(),0.0);
            volumeFractionBF_[i][j].setSize(patch.size(),0.0);
            numberFluxBF_[i][j].setSize(patch.size(),0.0);
            massFluxBF_[i][j].setSize(patch.size(),0.0);
            conductiveHeatFluxBF_[i][j].setSize(patch.size(),0.0);

        }
    }
    
    
    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidBoundaryMeasurements::clean()
{
    //- clean geometric fields
    
    forAll(rhoNBF_, i)
    {
        forAll(rhoNBF_[i], j)
        {
            
            momentumBF_[i][j] = Zero;
            UMeanBF_[i][j] = Zero;
            rhoNBF_[i][j] = 0.0;
            rhoMBF_[i][j] = 0.0;
            linearKEBF_[i][j] = 0.0;
            volumeFractionBF_[i][j] = 0.0;
            numberFluxBF_[i][j] = 0.0;
            massFluxBF_[i][j] = 0.0;
            conductiveHeatFluxBF_[i][j] = 0.0;
        }
    }
}


void solidBoundaryMeasurements::updateFields
(
    solidParticleCoupling& p
)
{
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
