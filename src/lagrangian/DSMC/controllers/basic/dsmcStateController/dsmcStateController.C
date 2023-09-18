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

#include "dsmcStateController.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dsmcStateController, 0);
defineRunTimeSelectionTable(dsmcStateController, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcStateController::dsmcStateController
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcControllerBase(mesh, cloud, dict),
    timePeriod_(timeDict_.get<scalar>("initialTimePeriod")),
    initialTime_(mesh_.time().startTime().value())
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find cellZone: " << regionName_ << nl << " in: "
            << mesh_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    const scalar avTimeInterval = timeData_.averageTimeInterval().deltaT();

    if ((timePeriod_ < avTimeInterval) && (timePeriod_ > 0.0))
    {
        timePeriod_ = avTimeInterval;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dsmcStateController> Foam::dsmcStateController::New
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("stateControllerModel"));

    Info<< "Selecting dsmcStateController " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "dsmcStateController",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<dsmcStateController>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcStateController::updateTime()
{
    dsmcControllerBase::updateTime();

    const scalar t = mesh_.time().timeOutputValue();

    if ((t - initialTime_) < timePeriod_)
    {
        timeData_.controlTimeInterval().endTime() = false;
    }
}


void Foam::dsmcStateController::updateProperties(const dictionary& dict)
{
    dsmcControllerBase::updateProperties(dict);

    timeDict_ = controllerDict_.subDict("timeProperties");
    timeDict_.readIfPresent("resetAtOutput", timeData_.resetFieldsAtOutput());
}


const Foam::labelList& Foam::dsmcStateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}


// ************************************************************************* //
