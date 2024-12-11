/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "diffusionBasedSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::diffusionBasedSmoothing::diffusionBasedSmoothing
(
    const dictionary& dict,
    const Time& runTime,
    const fvMesh& mesh,
    solidParcelCloud& cloud
)
:
    propertiesDict_(dict.subDict("diffusionBasedSmoothingProperties")),
    diffusionRunTime_
    (
        "diffisionControlDict",
        runTime.rootPath(),
        runTime.caseName()
    ),
    diffusionMesh_
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            diffusionRunTime_.timeName(),
            diffusionRunTime_,
            Foam::IOobject::MUST_READ
        )
    ),
    mesh_(mesh),
    simple_(diffusionMesh_),
    diffusionBandWidth_(propertiesDict_, "diffusionBandWidth", 0.01),
    diffusionSteps_(propertiesDict_, "diffusionSteps", 6),
    implicitFvm_(propertiesDict_, "useImplicitLaplacian", true),
    smoothDirection_
    (
        propertiesDict_.lookupOrDefault
        (
            "smoothDirection",
            tensor(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)
        )
    ),
    DT_("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_),
    startTime(diffusionRunTime_.startTime()),
    startTimeIndex(diffusionRunTime_.startTimeIndex()),
    diffusionTime_(0),
    diffusionDeltaT_(0)
    // diffusionBandWidth_(propertiesDict_.lookupOrDefault<scalar>("diffusionBandWidth", 0.006)),
    // diffusionSteps_(propertiesDict_.lookupOrDefault<scalar>("diffusionSteps", 6))
    
{

    // determine the time and time step in diffusion procedure
    diffusionTime_ = pow(diffusionBandWidth_.value(), 2)/4;
    diffusionDeltaT_ = diffusionTime_/diffusionSteps_.value();

    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionDeltaT_);
    Info << "Diffusion-based smoothing bandwidth is: "<< diffusionBandWidth_.value()<< nl
         << "Diffusion-based smoothing time is: " << diffusionRunTime_.endTime() << nl
         << "Diffusion-based smoothing time step is: " << diffusionRunTime_.deltaTValue() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::diffusionBasedSmoothing::~diffusionBasedSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::diffusionBasedSmoothing::diffusion
(
    volScalarField& s
)
{
    volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            dimless,
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
    );

    diffWorkField.ref().field() = s.internalField();

    Info<< "- Smoothing " << s.name() <<" -"<< endl;

    if(implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
            }
        }
    }
    else
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
            }
        }
    }
    
    diffusionRunTime_.setTime(startTime,startTimeIndex);

    s.ref().field() = diffWorkField.internalField();

    return;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::diffusionBasedSmoothing::diffusion
(
    const volScalarField::Internal& s
)
{    
    tmp<volScalarField::Internal> tS
    (
        new volScalarField::Internal
        (
            // "tS",
            IOobject
            (
                "tS",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", s.dimensions(), 0.0)
        )
    );

    scalarField& S = tS.ref();
    
    S = s;
    
    volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            s.dimensions(),
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
        
    );

    
    scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    diffWorkFieldInterFeildRef = S;

    Info<< "- Smoothing " << s.name() <<" -"<< endl;
    
    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
            }
        }

    }
    else
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
            }
        }
    }
    
    S = diffWorkField.internalField();
    
    diffusionRunTime_.setTime(startTime, startTimeIndex);

    return tS;
}

Foam::tmp<Foam::volVectorField::Internal> Foam::diffusionBasedSmoothing::diffusion
(
    const volVectorField::Internal& s
)
{ 
    tmp<volVectorField::Internal> tS
    (
        new volVectorField::Internal
        (
            // "tS",
            IOobject
            (
                "tS",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", s.dimensions(), vector::zero)
        )
    );

    vectorField& S = tS.ref();
    
    S = s;
    
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedVector
        (
            "zero",
            s.dimensions(),
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
        
    );

    
    vectorField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    diffWorkFieldInterFeildRef = S;

    Info<< "- Smoothing " << s.name() <<" -"<< endl;
    
    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
            }
        }

    }
    else
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
            }
        }
    }
    
    S = diffWorkField.internalField();
    
    diffusionRunTime_.setTime(startTime, startTimeIndex);

    return tS;
}

// ************************************************************************* //
