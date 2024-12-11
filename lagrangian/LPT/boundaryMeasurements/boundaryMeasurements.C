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
    boundaryMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "boundaryMeasurements.H"
#include "solidParcelCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from mesh and cloud 
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    solidParcelCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


//- Construct from mesh, cloud and boolean (dsmcFoam)
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    solidParcelCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIds_(identity(cloud_.typeIdList().size())),
    momentumBF_(),
    UMeanBF_(),
    fDBF_(),
    rhoNBF_(),
    rhoMBF_(),
    qBF_()
{
    
    momentumBF_.setSize(typeIds_.size());
    UMeanBF_.setSize(typeIds_.size());
    fDBF_.setSize(typeIds_.size());
    rhoNBF_.setSize(typeIds_.size());
    rhoMBF_.setSize(typeIds_.size());
    qBF_.setSize(typeIds_.size());
    
    forAll(rhoNBF_, i)
    {        
        momentumBF_[i].setSize(mesh_.boundaryMesh().size());
        UMeanBF_[i].setSize(mesh_.boundaryMesh().size());
        fDBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoMBF_[i].setSize(mesh_.boundaryMesh().size());
        qBF_[i].setSize(mesh_.boundaryMesh().size());
        
        forAll(rhoNBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            momentumBF_[i][j].setSize(patch.size(),vector::zero);
            UMeanBF_[i][j].setSize(patch.size(),vector::zero);
            fDBF_[i][j].setSize(patch.size(),vector::zero);
            rhoNBF_[i][j].setSize(patch.size(),0.0);
            rhoMBF_[i][j].setSize(patch.size(),0.0);
            qBF_[i][j].setSize(patch.size(),0.0);

        }
    }
    
    
    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void boundaryMeasurements::clean()
{
    //- clean geometric fields
    
    forAll(rhoNBF_, i)
    {
        forAll(rhoNBF_[i], j)
        {
            
            momentumBF_[i][j] = Zero;
            UMeanBF_[i][j] = Zero;
            fDBF_[i][j] = Zero;
            rhoNBF_[i][j] = 0.0;
            rhoMBF_[i][j] = 0.0;
            qBF_[i][j] = 0.0;
            
        }
    }
}


void boundaryMeasurements::updateFields
(
    solidParcel& p
)
{
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
