/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::diffusionBasedSmoothing

Group
    grpLagrangianIntermediateHeatTransferSubModels

Description
    Dummy heat transfer model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef diffusionBasedSmoothing_H
#define diffusionBasedSmoothing_H

#include "fvCFD.H"
#include "volMesh.H"
#include "simpleControl.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class solidParcelCloud;

/*---------------------------------------------------------------------------*\
                       Class diffusionBasedSmoothing Declaration
\*---------------------------------------------------------------------------*/
class diffusionBasedSmoothing
{
private:

    const dictionary propertiesDict_;

    Time diffusionRunTime_;

    fvMesh diffusionMesh_;

    const fvMesh& mesh_;

    simpleControl simple_;

     //- Phase diffusion band width
    demandDrivenEntry<scalar> diffusionBandWidth_;

    //- Diffusion number of time steps
    demandDrivenEntry<label> diffusionSteps_;

    //- Scheme True if use implicit laplacian, False if use explicit laplacian
    demandDrivenEntry<bool> implicitFvm_;

    tensor smoothDirection_;

    //- Diffusion tensor
    dimensionedTensor DT_;

    //- Diffusion start time
    dimensionedScalar startTime;
    label startTimeIndex;

    scalar diffusionTime_;

    scalar diffusionDeltaT_;

    
public:

    //- Runtime type information
    // TypeName("diffusionBasedSmoothing");


    // Constructors

        //- Construct from dictionary
        diffusionBasedSmoothing
        (
            const dictionary&, 
            const Time&, 
            const fvMesh& mesh,
            solidParcelCloud& cloud
        );

        //- Construct copy
        // diffusionBasedSmoothing(const diffusionBasedSmoothing& im);

        //- Construct and return a clone
        // virtual autoPtr<smoothingScheme> clone() const
        // {
        //     return autoPtr<smoothingScheme>
        //     (
        //         new diffusionBasedSmoothing(*this)
        //     );
        // }


    //- Destructor
    virtual ~diffusionBasedSmoothing();


    // Member Functions

        //- Return the diffused field obtained from the given volScalarField s
        void diffusion(volScalarField& s);

        //- Return the diffused field obtained from the given volVectorField s
        // void diffusion(volVectorField& s);

        //- Return the diffused field obtained from the given volScalarField s
        tmp<volScalarField::Internal> diffusion(const volScalarField::Internal& s);

        //- Return the diffused field obtained from the given volVectorField s
        tmp<volVectorField::Internal> diffusion(const volVectorField::Internal& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
