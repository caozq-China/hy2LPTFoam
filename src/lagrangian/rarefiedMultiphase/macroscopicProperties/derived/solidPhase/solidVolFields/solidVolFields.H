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

Class
    solidVolFields

Description
    Geometric field for temperature, density, velocity, pressure.



SourceFiles
    solidVolFieldsI.H
    solidVolFields.C
    solidVolFieldsIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidVolFields_H
#define solidVolFields_H

#include "solidField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidVolFields Declaration
\*---------------------------------------------------------------------------*/

class solidVolFields
:
    public solidField
{

private:

    // Private data
        dictionary propsDict_;
       
        //- Wall unit vectors - temporary storage, should be done in 
        // boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

        label sampleInterval_;
        label sampleCounter_;
//         scalar mfpReferenceTemperature_;
        word fieldName_;
        volScalarField solidRhoN_;
        volScalarField solidRhoNMean_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField loadBalanceWeight_;
        volScalarField overallTsolid_;
        volScalarField arithmeticAveT_;
        volScalarField surfaceNumberFlux_;
        volScalarField surfaceMassFlux_;
        volScalarField surfaceConductiveHeatFlux_;
        volScalarField liquidMassFraction_;
//         volScalarField volFraction_;
        volVectorField UMean_;
        label nTimeSteps_;
        labelList typeIdSolids_;
//         scalarField nTotColls_;
        scalarField rhoNMean_;
        scalarField rhoNInstantaneous_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMean_;
        scalarField rhoMMeanXnParticle_;
        scalarField liquidMassXnParticle_;
//         scalarField volumeXnParticle_;
        scalarField cXm_;
        scalarField cXmXT_;
        scalarField TXnParticle_;
        scalarField nColls_;
        vectorField momentumMean_;
        vectorField momentumMeanXnParticle_;
        IOList<scalar> nTotColls_;
        List<labelList> boundaryCells_;
        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        
        // boundary measurements
        
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
//         List<scalarField> volFractionBF_;
        List<scalarField> numberFluxBF_;
        List<scalarField> massFluxBF_;
        List<scalarField> conductiveHeatFluxBF_;
        List<vectorField> momentumBF_;
        
        bool averagingAcrossManyRuns_;
        bool densityOnly_;
//         bool measureHeatFluxShearStress_;
        
    // Private Member Functions

        void readIn();

        void writeOut();
       
        void calculateWallUnitVectors();


public:

    //- Runtime type information
    TypeName("solidVolFields");

    // Constructors

        //- Construct from components
    solidVolFields
    (
        const Time& t,
        const polyMesh& mesh,
        solidParticleCouplingCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    virtual ~solidVolFields() = default;


    // Member Functions

    //- create initial field
    virtual void createField();

    //- calculate field
    virtual void calculateField();

    //- write field
    virtual void writeField();
    
//     scalarField& rhoNInstantaneous();
//     volVectorField& UMean();
//     
//     volScalarField& rhoM();
    
//     word& fieldName();

//         const propertyField& fields() const;
    virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

