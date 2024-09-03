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
    particleCartesian

Description

\*----------------------------------------------------------------------------*/

#include "particleCartesian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleCartesian, 0);

    addToRunTimeSelectionTable
    (
        coordinateSystemType, 
        particleCartesian,
        fvMesh
    );
    
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
particleCartesian::particleCartesian
(
    Time& t,
    const polyMesh& mesh,
    solidParcelCloud& cloud
)
:
    coordinateSystemType(t, mesh, cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleCartesian::~particleCartesian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleCartesian::checkCoordinateSystemInputs(const bool init)
{
    writeCoordinateSystemInfo();
}


void particleCartesian::evolve()
{}


void particleCartesian::writeCoordinateSystemInfo() const
{}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
