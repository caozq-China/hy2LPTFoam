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

#include "solidStateController.H"
// #include "IFstream.H"
// #include "graph.H"
#include "solidParticleCouplingCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidStateController, 0);

defineRunTimeSelectionTable(solidStateController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidStateController::solidStateController
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
:
    solidControllerBase(mesh, cloud, dict),
    timePeriod_(timeDict_.get<scalar>("initialTimePeriod")),
    initialTime_(mesh_.time().startTime().value())
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << regionName_ << nl << "in: "
            << mesh_.time().system()/"solidControllersDict"
            << exit(FatalError);
    }
    
    const scalar avTimeInterval = timeData_.averageTimeInterval().deltaT();
    
    if ((timePeriod_ < avTimeInterval) && (timePeriod_ > 0.0))
    {
        timePeriod_ = avTimeInterval;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidStateController> solidStateController::New
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
{
    const word modelType
    (
        dict.lookup("stateControllerModel")
    );

    Info<< "Selecting solidStateController "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidStateController",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<solidStateController>
	(
		cstrIter()(mesh, cloud, dict)
	);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidStateController::updateTime()
{
    solidControllerBase::updateTime();

    const scalar t = mesh_.time().timeOutputValue();

    if ((t - initialTime_) < timePeriod_)
    {
        timeData_.controlTimeInterval().endTime() = false;
    }
}

void Foam::solidStateController::updateProperties(const dictionary& dict)
{
    solidControllerBase::updateProperties(dict);

    timeDict_ = controllerDict_.subDict("timeProperties");
    timeDict_.readIfPresent("resetAtOutput", timeData_.resetFieldsAtOutput());
}


const Foam::labelList& Foam::solidStateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}

} // End namespace Foam

// ************************************************************************* //
