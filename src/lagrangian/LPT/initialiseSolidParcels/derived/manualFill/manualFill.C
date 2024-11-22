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

#include "manualFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(manualFill, 0);

addToRunTimeSelectionTable(configuration, manualFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
manualFill::manualFill
(
    solidParcelCloud& cloud,
    const dictionary& dict
)
:
    configuration(cloud, dict),
    positionsFile_(dict.lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            cloud.time().constant(),
            cloud.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCells_(positions_.size(), -1),
    injectorTetFaces_(positions_.size(), -1),
    injectorTetPts_(positions_.size(), -1)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void manualFill::setInitialConfiguration()
{

    Info<< nl << "Initialising particles" << endl;

    const vector U(initialiseDict_.lookup("velocity"));

    forAll(positions_, i)
    {

        cloud_.mesh().findCellFacePt
        (
            positions_[i],
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i]
        );

        label typeId = 0;
        //const scalar& RWF = cloud_.coordSystem().RWF(injectorCells_[i]);
        
        scalar RWF = 1.0;
                        
        if(cloud_.axisymmetric())
        {                      
            const point& cC = cloud_.mesh().cellCentres()[injectorCells_[i]];
            scalar radius = cC.y();
            
            RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
        }

        cloud_.addNewParcel
        (
            mesh_,
            cloud_.constProps(typeId),
            positions_[i],
            U,
            RWF,
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i],
            typeId,
            -1
        );
    }
}


} // End namespace Foam

// ************************************************************************* //
