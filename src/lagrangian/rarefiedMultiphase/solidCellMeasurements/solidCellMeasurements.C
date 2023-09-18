/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "solidCellMeasurements.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidCellMeasurements::solidCellMeasurements
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    spc_(cloud)
{}


Foam::solidCellMeasurements::solidCellMeasurements
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const bool dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    spc_(cloud),
//     collisionSeparation_(),
    nColls_(),
    nCollsTotal_()
{
//     collisionSeparation_.setSize(mesh.nCells());
    nColls_.setSize(mesh.nCells());
    nCollsTotal_.setSize(Pstream::nProcs());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidCellMeasurements::clean()
{
    // Clean geometric fields
//     collisionSeparation_ = 0.0;
    nColls_ = 0.0;
    nCollsTotal_ = 0.0;
}


void Foam::solidCellMeasurements::updateFields(solidParticleCoupling& p)
{}


// ************************************************************************* //
