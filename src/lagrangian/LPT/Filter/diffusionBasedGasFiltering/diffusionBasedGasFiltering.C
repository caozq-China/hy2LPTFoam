/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "diffusionBasedGasFiltering.H"
#include "fvcLaplacianPlus.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::diffusionBasedGasFiltering::diffusionBasedGasFiltering
(
    const dictionary& dict, 
    const Time& runTime, 
    const fvMesh& mesh,
    solidParcelCloud& cloud
)
:
    propertiesDict_(dict.subDict("diffusionBasedGasFilteringProperties")),
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
    implicitFvm_(propertiesDict_, "useImplicitLaplacian", true),
    smoothDirection_
    (
        propertiesDict_.lookupOrDefault
        (
            "smoothDirection",
            tensor(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)
        )
    ),
    DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_),
    startTime(diffusionRunTime_.startTime()),
    startTimeIndex(diffusionRunTime_.startTimeIndex()),
    diffusionBandWidth_(propertiesDict_.lookupOrDefault<scalar>("diffusionBandWidth", 0.003)),
    diffusionSteps_(propertiesDict_.lookupOrDefault("diffusionSteps", 6)),
    diffusionTime_(0),
    diffusionDeltaT_(0),
    adjustDeltaT_(propertiesDict_, "adjustDiffusionSteps", true),
    stepFlux_(0),
    stepScalarField_(0),
    deltaTList_(0)
{    
    // determine the time and time step in diffusion procedure
    diffusionTime_ = pow(diffusionBandWidth_, 2)/4;
    diffusionDeltaT_ = diffusionTime_/diffusionSteps_;

    Info << "Diffusion-based gas filtering bandwidth is: "<< diffusionBandWidth_<< nl
         << "Diffusion-based gas filtering time is: " << diffusionTime_ << nl
         << "Diffusion-based gas filtering time step is: " << diffusionDeltaT_ << endl;

}
    
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// const access
inline Foam::scalar
Foam::diffusionBasedGasFiltering::diffusionBandWidth() const
{
    return diffusionBandWidth_;
}


inline Foam::label
Foam::diffusionBasedGasFiltering::diffusionSteps() const
{
    return diffusionSteps_;
}


inline Foam::scalar
Foam::diffusionBasedGasFiltering::diffusionDeltaT() const
{
    return diffusionDeltaT_;
}


inline Foam::scalar
Foam::diffusionBasedGasFiltering::diffusionTime() const
{
    return diffusionTime_;
}


// access to diffusion settings
inline Foam::scalar&
Foam::diffusionBasedGasFiltering::diffusionBandWidth()
{
    return diffusionBandWidth_;
}


inline Foam::label&
Foam::diffusionBasedGasFiltering::diffusionSteps()
{
    return diffusionSteps_;
}


inline Foam::scalar&
Foam::diffusionBasedGasFiltering::diffusionDeltaT()
{
    return diffusionDeltaT_;
}


inline Foam::scalar&
Foam::diffusionBasedGasFiltering::diffusionTime()
{
    return diffusionTime_;
}


//- update record list from explicit diffusion that need to reverse diffusion calculation      
void Foam::diffusionBasedGasFiltering::preExplicit(const volScalarField& s)
{
    // clear record
    stepFlux_.clear();
    stepScalarField_.clear();
    deltaTList_.clear();
    
    scalar diffusionExplicitDeltaT = diffusionDeltaT_;
    if (adjustDeltaT_.value())
    {
        diffusionExplicitDeltaT = diffusionDeltaT_/20;
    }
    label stepIndex = 1;// 1-10 diffusionDeltaT_/10; 11-20, diffusionDeltaT_/5, 
                        //>21 diffusionDeltaT_/2 >51 diffusionDeltaT_
    
    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionExplicitDeltaT);
    
    volScalarField S
    (
        IOobject
        (
            "tempExplicitDiffScalar",
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
    
    dimensionedTensor DTr("DTr", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);
    
    scalarField& diffWorkFieldInterFeildRef = S.ref();
    
    const scalarField& sInterFeildRef = s;

    diffWorkFieldInterFeildRef = sInterFeildRef;
    
    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            
        if (isA<processorPolyPatch>(pp))
        {
            scalarField& SBc = S.boundaryFieldRef()[patchi];
            SBc = s.boundaryField()[patchi];            
        }
    }
                
    while (diffusionRunTime_.loop())
    {
        while (simple_.correctNonOrthogonal())
        {
            
            fvScalarMatrix SEqn
            (
                fvm::ddt(S) 
             ==
                fvc::laplacian(DT, S)
            );
            
            stepFlux_.append(fvc::laplacianScheme2Grad(DTr, S, "gasFilterDiffusion2"));
            deltaTList_.append(diffusionRunTime_.deltaTValue());
          
            SEqn.solve();

            stepScalarField_.append(S);
        }
    
        stepIndex++;
        
        if (adjustDeltaT_.value())
        {
            scalar adjustDiffusionDeltaT = diffusionDeltaT_;
        
            if (stepIndex <= 4 )
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/20;
            }
            else if (stepIndex <= 10)
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/7.5;
            }
            else if (stepIndex <= 15)
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/5;
            }
            else if (stepIndex <= 21)
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/2;
            }
            else
            {
                adjustDiffusionDeltaT = diffusionDeltaT_;
            }
            
            diffusionRunTime_.setDeltaT(adjustDiffusionDeltaT);
        }
    }

    diffusionRunTime_.setTime(startTime, startTimeIndex);
}


void Foam::diffusionBasedGasFiltering::preExplicit
(
    const volScalarField& alpha,
    const volScalarField& s
)
{
    // clear record
    stepFlux_.clear();
    stepScalarField_.clear();
    deltaTList_.clear();
    
    scalar diffusionExplicitDeltaT = diffusionDeltaT_;
    if (adjustDeltaT_.value())
    {
        diffusionExplicitDeltaT = diffusionDeltaT_/20;
    }
    label stepIndex = 1;// 1-10 diffusionDeltaT_/10; 11-20, diffusionDeltaT_/5, 
                        //>21 diffusionDeltaT_/2 >51 diffusionDeltaT_
    
    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionExplicitDeltaT);
    
    volScalarField S
    (
        IOobject
        (
            "tempExplicitDiffScalar",
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
    
    dimensionedTensor DTr("DTr", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);
    
    scalarField& diffWorkFieldInterFeildRef = S.ref();
    
    const scalarField& sInterFeildRef = s;

    diffWorkFieldInterFeildRef = sInterFeildRef;
    
    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            
        if (isA<processorPolyPatch>(pp))
        {
            scalarField& SBc = S.boundaryFieldRef()[patchi];
            SBc = s.boundaryField()[patchi];
        }
    }
       
    
    while (diffusionRunTime_.loop())
    {
        while (simple_.correctNonOrthogonal())
        {
            fvScalarMatrix SEqn
            (
                fvm::ddt(S) 
             ==
                fvc::laplacian(fvc::interpolate(alpha)*DT, S)
            );
            
            stepFlux_.append(fvc::laplacianScheme2Grad(fvc::interpolate(alpha)*DTr, S, "gasFilterDiffusion2"));
            deltaTList_.append(diffusionRunTime_.deltaTValue());

            SEqn.solve();
    
            stepScalarField_.append(S);
        }
        
        stepIndex++;
        
        if (adjustDeltaT_.value())
        {
            scalar adjustDiffusionDeltaT = diffusionDeltaT_;
        
            if (stepIndex <= 4 )
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/20;
            }
            else if (stepIndex <= 10)
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/7.5;
            }
            else if (stepIndex <= 15)
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/5;
            }
            else if (stepIndex <= 21)
            {
                adjustDiffusionDeltaT = diffusionDeltaT_/2;
            }
            else
            {
                adjustDiffusionDeltaT = diffusionDeltaT_;
            }
            
            diffusionRunTime_.setDeltaT(adjustDiffusionDeltaT);
        }
    }
    
    diffusionRunTime_.setTime(startTime, startTimeIndex);
}



//- reverse explicit laplacian diffusion calculation
tmp<volScalarField>
Foam::diffusionBasedGasFiltering::reverseExplicit(const volScalarField& F)
{
    const label stopTimeIndex = deltaTList_.size() - 1;
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const scalarField& meshVolume = mesh_.V();
 
    tmp<volScalarField> filteredField
    (
        new volScalarField
        (
            IOobject
            (
                "gasFilteredField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", F.dimensions(), 0.0)
        )
    );

    volScalarField& field = filteredField.ref();
    field = F;
    field = field*stepScalarField_[stopTimeIndex];
    scalarField& iSS = field.primitiveFieldRef();

    label indexReverse;
    
    forAll (stepFlux_, timeIndex)
    {
        indexReverse = stopTimeIndex - timeIndex;
        
        const scalarField iSSref = field;

        scalar deltaT = deltaTList_[indexReverse];
    
        forAll(owner, facei)
        {
            
            if (stepScalarField_[indexReverse][owner[facei]] > SMALL && stepScalarField_[indexReverse][neighbour[facei]] > SMALL)
            {
                if (stepFlux_[indexReverse][facei]>0)
                {
                    iSS[owner[facei]] -= iSSref[owner[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[owner[facei]]/stepScalarField_[indexReverse][owner[facei]];
                    iSS[neighbour[facei]] += iSSref[owner[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[owner[facei]]/stepScalarField_[indexReverse][owner[facei]];
                }
                else
                {
                    iSS[owner[facei]] -= iSSref[neighbour[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[neighbour[facei]]/stepScalarField_[indexReverse][neighbour[facei]];
                    iSS[neighbour[facei]] += iSSref[neighbour[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[neighbour[facei]]/stepScalarField_[indexReverse][neighbour[facei]];
                }
                
            }
            
        }      
        
        forAll(mesh_.boundary(), patchi)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            
            if (isA<processorPolyPatch>(pp))
            {
                
                const labelUList& pFaceCells = mesh_.boundary()[patchi].faceCells();

                const fvsPatchScalarField& pssf = stepFlux_[indexReverse].boundaryField()[patchi];

                forAll(mesh_.boundary()[patchi], facei)
                {
                    if (stepScalarField_[indexReverse][pFaceCells[facei]] > SMALL && stepScalarField_[indexReverse].boundaryField()[patchi][facei] > SMALL)
                    {                        
                        if (pssf[facei]>0)
                        {
                            iSS[pFaceCells[facei]] -= iSSref[pFaceCells[facei]]*pssf[facei]*deltaT/meshVolume[pFaceCells[facei]]/stepScalarField_[indexReverse][pFaceCells[facei]];
                        }
                        else
                        {
                            iSS[pFaceCells[facei]] -= field.boundaryField()[patchi][facei]*pssf[facei]*deltaT/meshVolume[pFaceCells[facei]]/stepScalarField_[indexReverse].boundaryField()[patchi][facei];
                        }
                    }
                }
            }
        }
        
        field.correctBoundaryConditions();
    }        
    return filteredField;
}


tmp<volVectorField>
Foam::diffusionBasedGasFiltering::reverseExplicit(const volVectorField& F)
{
    const label stopTimeIndex = deltaTList_.size() - 1;
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const scalarField& meshVolume = mesh_.V();
    
    tmp<volVectorField> filteredField
    (
        new volVectorField
        (
            IOobject
            (
                "gasFilteredField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", F.dimensions(), vector::zero)
        )
    );

    volVectorField& field = filteredField.ref();
    field = F;
    field = field*stepScalarField_[stopTimeIndex];
    vectorField& iSS = field.primitiveFieldRef();
    
    label indexReverse;

    forAll (stepFlux_, timeIndex)
    {
        indexReverse = stopTimeIndex - timeIndex;
        
        const vectorField iSSref = field;
        
        scalar deltaT = deltaTList_[indexReverse];
        
    
        forAll(owner, facei)
        {
            
            if (stepScalarField_[indexReverse][owner[facei]] > SMALL && stepScalarField_[indexReverse][neighbour[facei]] > SMALL)
            {
                if (stepFlux_[indexReverse][facei]>0)
                {
                    iSS[owner[facei]] -= iSSref[owner[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[owner[facei]]/stepScalarField_[indexReverse][owner[facei]];
                    iSS[neighbour[facei]] += iSSref[owner[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[owner[facei]]/stepScalarField_[indexReverse][owner[facei]];
                }
                else
                {
                    iSS[owner[facei]] -= iSSref[neighbour[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[neighbour[facei]]/stepScalarField_[indexReverse][neighbour[facei]];
                    iSS[neighbour[facei]] += iSSref[neighbour[facei]]*stepFlux_[indexReverse][facei]*deltaT/meshVolume[neighbour[facei]]/stepScalarField_[indexReverse][neighbour[facei]];
                }
                
            }
            
        }
        forAll(mesh_.boundary(), patchi)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            
            if (isA<processorPolyPatch>(pp))
            {
                
                const labelUList& pFaceCells = mesh_.boundary()[patchi].faceCells();

                const fvsPatchScalarField& pssf = stepFlux_[indexReverse].boundaryField()[patchi];

                forAll(mesh_.boundary()[patchi], facei)
                {
                    if (stepScalarField_[indexReverse][pFaceCells[facei]] > SMALL && stepScalarField_[indexReverse].boundaryField()[patchi][facei] > SMALL)
                    { 
                        if (pssf[facei]>0)
                        {
                            iSS[pFaceCells[facei]] -= iSSref[pFaceCells[facei]]*pssf[facei]*deltaT/meshVolume[pFaceCells[facei]]/stepScalarField_[indexReverse][pFaceCells[facei]];
                        }
                        else
                        {
                            iSS[pFaceCells[facei]] -= field.boundaryField()[patchi][facei]*pssf[facei]*deltaT/meshVolume[pFaceCells[facei]]/stepScalarField_[indexReverse].boundaryField()[patchi][facei];
                        }
                    }
                }
            }
        }
        
        field.correctBoundaryConditions();
    }
    
    return filteredField;
}



// Return the filtered field from the given volScalarField F
tmp<volScalarField>
Foam::diffusionBasedGasFiltering::filteredField(const volScalarField& F)
{
    if (implicitFvm_.value())
    {
        
        diffusionRunTime_.setEndTime(diffusionTime_);
        diffusionRunTime_.setDeltaT(diffusionDeltaT_);
        
        tmp<volScalarField> tF
        (
            new volScalarField
            (
                IOobject
                (
                    "tF",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", F.dimensions(), 0.0)
            )
        );

        volScalarField& S = tF.ref();
        
        S = F;
        
        scalarField& iS = S;
        
        volScalarField diffWorkField
        (
            IOobject
            (
                "tempGasDiffScalar",
                diffusionRunTime_.timeName(),
                diffusionMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            diffusionMesh_,
            dimensionedScalar
            (
                "zero",
                S.dimensions(),
                scalar(0.0)
            ),
            zeroGradientFvPatchScalarField::typeName
            
        );

        
        scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();
        
        diffWorkFieldInterFeildRef = iS;
        
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField,"gasFilterDiffusion"));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField,"gasFilterDiffusion"));
            }
        }
        
        iS = diffWorkFieldInterFeildRef;
        
        S.correctBoundaryConditions();
        
        diffusionRunTime_.setTime(startTime, startTimeIndex);
        
        return tF;

    }
    else
    {
        return reverseExplicit(F);
    }
}


tmp<volVectorField>
Foam::diffusionBasedGasFiltering::filteredField(const volVectorField& F)
{
    if (implicitFvm_.value())
    {
        
        diffusionRunTime_.setEndTime(diffusionTime_);
        diffusionRunTime_.setDeltaT(diffusionDeltaT_);
        
        tmp<volVectorField> tF
        (
            new volVectorField
            (
                // "tF",
                IOobject
                (
                    "tF",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", F.dimensions(), vector::zero)
            )
        );

        volVectorField& S = tF.ref();
        
        S = F;
        
        vectorField& iS = S;
        
        volVectorField diffWorkField
        (
            IOobject
            (
                "tempGasDiffVector",
                diffusionRunTime_.timeName(),
                diffusionMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            diffusionMesh_,
            dimensionedVector
            (
                "zero",
                S.dimensions(),
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
                
        );

        
        vectorField& diffWorkFieldInterFeildRef = diffWorkField.ref();
        
        diffWorkFieldInterFeildRef = iS;
        
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField,"gasFilterDiffusion"));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField,"gasFilterDiffusion"));
            }
        }
        
        iS = diffWorkFieldInterFeildRef;
        
        S.correctBoundaryConditions();
        
        diffusionRunTime_.setTime(startTime, startTimeIndex);
        
        return tF;

    }
    else
    {
        return reverseExplicit(F);
    }
}


tmp<volScalarField>
Foam::diffusionBasedGasFiltering::filteredField
(
    const volScalarField& alpha,
    const volScalarField& F
)
{
    if (implicitFvm_.value())
    {
        
        diffusionRunTime_.setEndTime(diffusionTime_);
        diffusionRunTime_.setDeltaT(diffusionDeltaT_);
        
        tmp<volScalarField> tF
        (
            new volScalarField
            (
                IOobject
                (
                    "tF",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", F.dimensions(), 0.0)
            )
        );

        volScalarField& S = tF.ref();
        
        S = F;
        
        scalarField& iS = S;
        
        volScalarField diffWorkField
        (
            IOobject
            (
                "tempGasDiffScalar",
                diffusionRunTime_.timeName(),
                diffusionMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            diffusionMesh_,
            dimensionedScalar
            (
                "zero",
                S.dimensions(),
                scalar(0.0)
            ),
            zeroGradientFvPatchScalarField::typeName
            
        );

        
        scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();
        
        diffWorkFieldInterFeildRef = iS;
        
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(fvc::interpolate(alpha)*DT, diffWorkField,"gasFilterDiffusion"));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(fvc::interpolate(alpha)*DT, diffWorkField,"gasFilterDiffusion"));
            }
        }
        
        iS = diffWorkFieldInterFeildRef;
        
        S.correctBoundaryConditions();
        
        diffusionRunTime_.setTime(startTime, startTimeIndex);
        
        return tF;

    }
    else
    {
        return reverseExplicit(F);
    }
}


tmp<volVectorField>
Foam::diffusionBasedGasFiltering::filteredField
(
    const volScalarField& alpha,
    const volVectorField& F
)
{
    if (implicitFvm_.value())
    {
        
        diffusionRunTime_.setEndTime(diffusionTime_);
        diffusionRunTime_.setDeltaT(diffusionDeltaT_);
        
        tmp<volVectorField> tF
        (
            new volVectorField
            (
                // "tF",
                IOobject
                (
                    "tF",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", F.dimensions(), vector::zero)
            )
        );

        volVectorField& S = tF.ref();
        
        S = F;
        
        vectorField& iS = S;
        
        volVectorField diffWorkField
        (
            IOobject
            (
                "tempGasDiffvector",
                diffusionRunTime_.timeName(),
                diffusionMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            diffusionMesh_,
            dimensionedVector
            (
                "zero",
                S.dimensions(),
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
                
        );

        
        vectorField& diffWorkFieldInterFeildRef = diffWorkField.ref();
        
        diffWorkFieldInterFeildRef = iS;
        
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(diffWorkField) - fvm::laplacian(fvc::interpolate(alpha)*DT, diffWorkField,"gasFilterDiffusion"));
                }
            }
            else
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(fvc::interpolate(alpha)*DT, diffWorkField,"gasFilterDiffusion"));
            }
        }
        
        iS = diffWorkFieldInterFeildRef;
        
        S.correctBoundaryConditions();
        
        diffusionRunTime_.setTime(startTime, startTimeIndex);
        
        return tF;

    }
    else
    {
        return reverseExplicit(F);
    }
}





// ************************************************************************* //
