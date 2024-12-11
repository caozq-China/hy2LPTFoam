/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

Class
    constantTimeStepModel

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "constantTimeStepModel.H"
#include "solidParcelCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantTimeStepModel, 0);

    addToRunTimeSelectionTable
    (
        TimeStepModel,
        constantTimeStepModel,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
constantTimeStepModel::constantTimeStepModel
(
    Time& t,
    const polyMesh& mesh,
    solidParcelCloud& cloud
)
:
    TimeStepModel(t, mesh, cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantTimeStepModel::~constantTimeStepModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void constantTimeStepModel::checkTimeStepModelInputs()
{
    writeTimeStepModelInfo();
}


void constantTimeStepModel::update()
{}


void constantTimeStepModel::writeTimeStepModelInfo() const
{
    Info<< "Constant time-step model:" << nl
        << "- time-step [sec]" << tab << deltaTValue() << nl
        << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
