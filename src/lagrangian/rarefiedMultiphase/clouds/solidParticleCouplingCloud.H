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

Class
    Foam::solidParticleCouplingCloud

Description
    A Cloud of solid particles

SourceFiles
    solidParticleCouplingCloudI.H
    solidParticleCouplingCloud.C
    solidParticleCouplingCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticleCouplingCloud_H
#define solidParticleCouplingCloud_H

#include "CloudWithModels.H"
//#include "Cloud.H"

#include "solidParticleCoupling.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"

#include "solidBoundaries.H"
#include "solidAllConfigurations.H"
#include "solidFieldProperties.H"

#include "dsmcCloud.H"
#include "InterphaseCoupling.H"
#include "solidPhaseChangeModel.H"

#include "solidBoundaryMeasurements.H"
#include "solidCellMeasurements.H"
#include "Tensor.H"

#include "solidCollisionDetection.H"
#include "solidBinaryCollisionModel.H"
#include "solidControllers.H"

// #include "surfaceErosionModel.H"
#include "mppicPackingModel.H"
#include "mppicDampingModel.H"
#include "mppicIsotropyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class mppicPackingModel;
class mppicDampingModel;
class mppicIsotropyModel;

/*---------------------------------------------------------------------------*\
                           Class solidParticleCouplingCloud Declaration
\*---------------------------------------------------------------------------*/

class solidParticleCouplingCloud
:
      public CloudWithModels<solidParticleCoupling>
{
public:

    typedef solidParticleCouplingCloud baseCloudType;
    
private:
        
    // Private data
        IOdictionary mppicPropertiesDict_;
        
        List<word> typeIdSolidList_;
        
        List<word> materialList_;
        
        //- number of actual solid particles represented by a computational solid particle
        scalar nSolidParticles_;
        
        scalar nDsmcReflectionCandidate_;
        
        scalar solidWeightFactor_;
        
        scalar criticalVolFrac_;
        
        bool enableMPPICMethod_;

//         bool enableSurfaceErosionModel_;
        
        volScalarField sigmaTcRMax_;
        
        //- for particle-particle collision model
        volScalarField cRsigmaPPMax_;
        
        volScalarField volumeFraction_;
        
//         volScalarField massBelowLunarSurface_;
        
//         volScalarField massErosionRate_;
        
        List<solidParticleCoupling::constantProperties> constSolidProps_;
        
        Random rndGenS_;
        
        //- controllers
        solidControllers controllers_;
        
        solidFieldProperties fields_; 
        
        //- boundaries
        solidBoundaries boundaries_;
        
        autoPtr<InterphaseCoupling> interphaseCouplingModel_;
        
        //- for selecting particle-particle collision partners
        autoPtr<solidCollisionDetection> solidCollisionDetectionModel_;
        
        //- for particle-particle collision
        autoPtr<solidBinaryCollisionModel> solidBinaryCollisionModel_;
        
        autoPtr<solidPhaseChangeModel> solidPhaseChange_;
        
//         autoPtr<surfaceErosionModel> surfaceErosionModel_;
        
        //- MPPIC
        DynamicList<label> cellSolidPackingOccupIds_;
        DynamicList<label> cellSolidCollisionOccupIds_;
        //- generate a instance of Packing model
        autoPtr<mppicDampingModel> mppicDampingModel_;
        autoPtr<mppicPackingModel> mppicPackingModel_;
        autoPtr<mppicIsotropyModel> mppicIsotropyModel_;

//         autoPtr<List<solidParticleCoupling*>> particlePtr_;
        
        solidBoundaryMeasurements solidBoundaryMeas_;
        
        solidCellMeasurements solidCellMeas_;

        dsmcCloud* dsmcCloudReference_;
        
        
        
    // Private Member Functions

        bool enableMoleculeVelocityUpdate();
        
//         bool enableStochasticCollision();
        
        void interphaseCouplingCal();
        
        void moleculePostCollisionVelocityUpdate
        (
            solidParticleCoupling& pP,
            dsmcParcel& pG
        );
        
        void buildSolidConstProps();
        
        void updateVelocityHalfTimeStep();
        
        void updateVelocityHalfTimeStepBack();
        
        void updateVelocityFullTimeStep();
        
        void updateMppicCellList();
        
        void clearMppicCellList();
        
        void particleParticleCollisions();
        
        //- Disallow default bitwise copy construct
        solidParticleCouplingCloud(const solidParticleCouplingCloud&);
        
        //- Disallow default bitwise assignment
        void operator=(const solidParticleCouplingCloud&);

public:
    
    // Constructors

        //- Construct given name and mesh, will read particles, DSMC parcels
        //- and fields from files
        solidParticleCouplingCloud
        (
            const Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            dsmcCloud* dsmcCloudReference,
            bool readFields = true
        );
        
        //- construct given name, mesh and initialisation dictionary.
        solidParticleCouplingCloud
        (
            const Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& solidInitialiseDict,
            dsmcCloud* dsmcCloudReference, 
            const bool& clearFields
        );
        
    virtual ~solidParticleCouplingCloud() = default;
    // Member Functions

        // Access
            labelList getTypeIDs(const dictionary& dict) const;
    
            inline const IOdictionary& mppicProperties() const;
            //- Return the particle IDs corresponding to the particle names read
            //- from the input dictionary
            
            inline const List<word>& typeIdSolidList() const;
            
            inline const List<word>& materialList() const;
            
            inline scalar nSolidParticles() const;
            
            inline scalar solidWeightingAmplificationFactor() const;
            
            inline scalar criticalVolFrac() const;
            
            inline bool enableMPPICMethod();
            
//             inline bool enableSurfaceErosionModel();
            
            inline volScalarField& sigmaTcRMax();
            
            inline volScalarField& cRsigmaPPMax();
            
            inline volScalarField& volumeFrac();
//             inline volScalarField& massBelowLunarSurface();
            
//             inline volScalarField& massErosionRate();
            
            //- Return all of the constant properties
            inline const List<solidParticleCoupling::constantProperties>& constSolidProps() const;

            //- Return the constant properties of the given typeId
            inline const solidParticleCoupling::constantProperties& constSolidProps(label typeIdSolid) const;
            
            //- Return reference to the random object 
            inline Random& rndGenS();
            
        //- Sub-models
            //- Return reference to binary elastic collision model
            inline const solidBinaryCollisionModel&  solidBinaryCollision() const;
            //- Return non-const reference to binary elastic collision model
            inline solidBinaryCollisionModel& solidBinaryCollision();
            
            inline const solidCollisionDetection&  solidCollisionDetectionModel() const;
            //- Return non-const reference to binary elastic collision model
            inline solidCollisionDetection& solidCollisionDetectionModel();
            
            
            inline const solidBoundaries& boundaries() const;  
            inline solidBoundaries& boundaries();
            
            inline const solidBoundaryMeasurements& boundaryFluxMeasurements() const;
            inline solidBoundaryMeasurements& boundaryFluxMeasurements();
            
            inline const solidCellMeasurements& solidCellPropMeasurements() const;
            inline solidCellMeasurements& solidCellPropMeasurements();
            
            inline const solidControllers& controllers() const;  //****
            inline solidControllers& controllers();//****
            
            //- interphase coupling models
            inline const InterphaseCoupling& interphaseCoupling() const;
            inline InterphaseCoupling& interphaseCoupling();
            
            //- phase change models
            inline const solidPhaseChangeModel& solidPhaseChange() const;
            inline solidPhaseChangeModel& solidPhaseChange(); 

//             inline const surfaceErosionModel& surfaceErosionModels() const;
//             inline surfaceErosionModel& surfaceErosionModels();
            inline const DynamicList<label>& cellSolidPackingOccupancyIds() const;
            inline const DynamicList<label>& cellSolidCollisionOccupancyIds() const;
            
            inline const mppicDampingModel& mppicDampingModels() const;
            inline mppicDampingModel& mppicDampingModels(); 
            
            inline const mppicPackingModel& mppicPackingModels() const;
            inline mppicPackingModel& mppicPackingModels(); 
            
            inline const mppicIsotropyModel& mppicIsotropyModels() const;
            inline mppicIsotropyModel& mppicIsotropyModels(); 
            
            //- Return the dsmcCloud 
            inline dsmcCloud* dsmcCloudReference() const;
            
//             inline label particlePhaseStateCheck
//             (
//                 const scalar temperature,
//                 const label typeID,
//                 const label CzRatio
//                 
//             );
            
            //- Particle diameter correction due to phase change
//             inline scalar particleDiameterCorrection
//             (
//                 const label phaseState,
//                 const label typeID,
//                 const label CzRatio
//             );
            
        
        //- check
            scalar conductiveEnergyTransferWithWallPatch
            (
                const solidParticleCoupling& p,
                const scalar Twall,
                const scalar poissonRatioWall,
                const scalar elasticModuliWall,
                const scalar materialDensityWall,
                const scalar CpWall,
                const scalar kWall
            );
            
            scalar particleLiquidDensityCorrection
            (
                const word materialName,
                const scalar temperature  
            ) const;
            
            scalar particleLatentHeatOfFusion
            (
                const word materialName
            ) const;
            
            scalar phaseChangeModelConst
            (
                const word materialName
            ) const;
            
//             scalar dParcel
//             (
//                 solidParticleCoupling& p
//             );

            void resetUcorrect();

            void updateParticleNumberSequence();
            
            void MPPICprocedures
            (
                solidParticleCoupling::trackingData& tdSolid
            );
            
            //- Cloud evolution functions
            void evolve();
            
            void axisymmetricWeighting();
            
//             void phaseCheck(solidParticleCoupling& pSolid);
            void particleSizeCorrection
            (
                solidParticleCoupling& pSolid
            );
            
            void updateParticleVolumeFraction();
            
            void addNewSolidParticle
            (
                const vector& position,
                const vector& U,
                const vector& omega,
                const vector& UCorrect,
                const vector& F,
                const scalar D,
                const scalar T,
                const scalar RWF,
                const scalar CzRatio,
                const label cellI,
                const label numSeq,
                const label typeId,
                const label phaseState,
                const label newSolidParticle
            );
            
            void info() const;
            
//             virtual bool hasWallImpactDistance() const;
            //- clear all existing solid particles in the field
            inline void clearSolidParticles();
            
//             CloudFunctionObjectList<baseCloudType>& functions()
//             {
//                 return functions_;
//             }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleCouplingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
