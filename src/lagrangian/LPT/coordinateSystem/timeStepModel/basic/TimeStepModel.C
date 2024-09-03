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
    TimeStepModel

Description

\*----------------------------------------------------------------------------*/

#include "TimeStepModel.H"
#include "solidParcelCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(TimeStepModel, 0);
    defineRunTimeSelectionTable(TimeStepModel, fvMesh);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::TimeStepModel::initialisenParticles(const scalar value)
{
    forAll(nParticles_, celli)
    {
        nParticles_[celli] = value;
    }

    forAll(nParticles_.boundaryField(), patchi)
    {
        forAll(nParticles_.boundaryField()[patchi], facei)
        {
            nParticles_.boundaryFieldRef()[patchi][facei] = value;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
Foam::TimeStepModel::TimeStepModel
(
    Time& t,
    const polyMesh& mesh,
    solidParcelCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    nParticles_
    (
        IOobject
        (
            "nParticles",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("nParticles", dimless, 1.0)
    )
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::TimeStepModel>
Foam::TimeStepModel::New
(
    Time& t,
    const polyMesh& mesh,
    solidParcelCloud& cloud
)
{
    word timeStepModel =cloud.particleProperties().lookup("timeStepModel");
        // cloud.particleProperties().lookupOrDefault<word>
        // (
        //     "timeStepModel",
        //     "constant"
        // );

    // timeStepModel = "dsmc" + static_cast<word>(std::toupper(timeStepModel[0]))
        // + timeStepModel.substr(1) + "TimeStepModel";

    Info<< "Selecting the time-step model:" << tab << timeStepModel
        << "\n" << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(timeStepModel);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "TimeStepModel::New(Time&, const polyMesh&, solidParcelCloud&)"
        )   << "Unknown time-step model type "
            << timeStepModel << endl << endl
            << "Valid time-step model types are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<TimeStepModel>(cstrIter()(t, mesh, cloud));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TimeStepModel::~TimeStepModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
