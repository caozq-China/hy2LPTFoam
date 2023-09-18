/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "solidGravitationalAccelerationController.H"
#include "addToRunTimeSelectionTable.H"
// #include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(solidGravitationalAccelerationController, 0);

addToRunTimeSelectionTable(solidStateController, solidGravitationalAccelerationController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidGravitationalAccelerationController::solidGravitationalAccelerationController
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
:
solidStateController(mesh, cloud, dict),
propsDict_(dict.subDict(typeName + "Properties")),
acceleration_(propsDict_.get<vector>("acceleration"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    singleValueController() = true;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidGravitationalAccelerationController::initialConfiguration()
{}

void solidGravitationalAccelerationController::calculateProperties()
{}

void solidGravitationalAccelerationController::controlSolidParticlesBeforeMove()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label cellI = controlZone()[c];

        for (solidParticleCoupling* p : cellOccupancy[cellI])
        {
            p->U() += 0.5*acceleration_*mesh_.time().deltaTValue();
        }
    }
}

void solidGravitationalAccelerationController::controlSolidParticlesAfterMove()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label cellI = controlZone()[c];

        for (solidParticleCoupling* p : cellOccupancy[cellI])
        {
            p->U() += 0.5*acceleration_*mesh_.time().deltaTValue();
        }
    }
}

void solidGravitationalAccelerationController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = mesh_.time();
    if(runTime.outputTime())
    {
//         NotImplemented;
    }
}

void solidGravitationalAccelerationController::controlSolidParticlesBeforeCollisions()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label cellI = controlZone()[c];

        for (solidParticleCoupling* p : cellOccupancy[cellI])
        {
            p->U() -= 0.5*acceleration_*mesh_.time().deltaTValue();
        }
    }
}

void solidGravitationalAccelerationController::controlSolidParticlesAfterCollisions()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label cellI = controlZone()[c];

        for (solidParticleCoupling* p : cellOccupancy[cellI])
        {
            p->U() += acceleration_*mesh_.time().deltaTValue();
        }
    }
}

void solidGravitationalAccelerationController::updateProperties(const dictionary& dict)
{
    solidStateController::updateProperties(dict);
    propsDict_ = dict.subDict(typeName + "Properties");
    setProperties();
}

void solidGravitationalAccelerationController::setProperties()
{
    acceleration_ = propsDict_.get<vector>("acceleration");
}


}

// ************************************************************************* //
