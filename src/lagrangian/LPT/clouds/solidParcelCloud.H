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
    Foam::solidParcelCloud

Description
    A Cloud of solid particles

SourceFiles
    solidParcelCloudI.H
    solidParcelCloud.C
    solidParcelCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParcelCloud_H
#define solidParcelCloud_H


#include "solidParcel.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "integrationScheme.H"
#include "cloudSolution.H"
#include "forceList.H"
#include "Tensor.H"
#include "solidBoundaries.H"
#include "allConfigurations.H"
#include "rho2ReactionThermo.H"
#include "HeatTransferModel.H"
#include "boundaryMeasurements.H"
#include "solidFieldProperties.H"
//#include "coordinateSystemType.H"
#include "diffusionBasedSmoothing.H"
#include "diffusionBasedGasFiltering.H"

#include "PackingModel.H"
#include "DampingModel.H"
#include "IsotropyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// class integrationScheme;

class PackingModel;
class DampingModel;
class IsotropyModel;

/*---------------------------------------------------------------------------*\
                           Class solidParcelCloud Declaration
\*---------------------------------------------------------------------------*/

class solidParcelCloud
:
      public Cloud<solidParcel>
{
private:
        
    // Private data

        const word& cloudName_;

        const fvMesh& mesh_;

        //- Dictionary of particle properties
        const IOdictionary particleProperties_;

        //- Dictionary of simulation properties
        // const dictionary controlDict_;

        //- Sub-models dictionary
        const dictionary subModelProperties_;

        List<word> typeIdList_;

        //- Parcel constant properties - one for each type
        List<solidParcel::constantProperties> constProps_;

        //- Gravity
        const dimensionedVector& g_;

        scalar nParticle_;
        scalar radialExtent_;
        scalar maxRWF_;

        label interiorInteractionType_;

        //- Include radiation
        // Switch radiation_;

        // Switch UreduceD_;
        bool axisymmetric_;
        Switch PreduceD_;
        Switch hasWallImpactDistance_;

        //- SLG thermodynamics package
        const rho2ReactionThermo& thermo_;

        // References to the gas fields
        const volVectorField& U_;//- gas Velocity [m/s]
        const volScalarField& Ma_;//- gas Mach number
        const volScalarField& gamma_;//- gas specific heat ratio
        const volScalarField& omega_;//- gas viscosity index

        volVectorField UFilter_;//- gas Velocity [m/s]
        volScalarField rhoFilter_;//- gas Density [kg/m3]
        volScalarField muFilter_;//- gas  Dynamic viscosity [Pa.s]
        volScalarField TtrFilter_;//- gas Temperature [K]
        volScalarField CpFilter_;//- gas specific translational-rotational heat capacity
        volScalarField kappaFilter_;//- gas heat conductivity
        volScalarField MaFilter_;//- gas Mach number
        volScalarField gammaFilter_;//- gas specific heat ratio
        volScalarField omegaFilter_;//- gas viscosity index

        volScalarField theta_;//- particle volume fraction
        volScalarField np_;//- number of real particles

        //- Cell length scale
        scalarField cellLengthScale_;

        //- Coordinate system model: Cartesian, axisymmetric, spherical
        //autoPtr<coordinateSystemType> particleCoordinateSystem_;
        
        //- A data structure holding which particles are in which cell
        autoPtr<List<DynamicList<solidParcel*>>> cellOccupancyPtr_;

    // models

        Random rndGen_;

        //- Solution properties
        cloudSolution solution_;

        forceList forcesList_;

        solidBoundaries boundaries_;

        boundaryMeasurements boundaryMeas_;

        solidFieldProperties fieldProps_;

        diffusionBasedSmoothing DBS_;
        diffusionBasedGasFiltering DBGF_;

        // Reference to the particle integration schemes
        //- Velocity integration
        autoPtr<integrationScheme> UIntegrator_;
        //- Temperature integration
        autoPtr<integrationScheme> TIntegrator_;

        //- MPPIC model
        autoPtr<PackingModel> packingModel_;
        autoPtr<DampingModel> dampingModel_;
        autoPtr<IsotropyModel> isotropyModel_;

        autoPtr<HeatTransferModel> heatTransferModel_;


    // Source

        //- Momentum
        autoPtr<volVectorField::Internal> momentumTrans_;

        //- Coefficient for carrier phase U equation
        autoPtr<volScalarField::Internal> momentumCoeff_;

        //- Sensible enthalpy transfer [J/kg]
        autoPtr<volScalarField::Internal> qTrans_;

        //- Coefficient for carrier phase hs equation [W/K]
        autoPtr<volScalarField::Internal> qCoeff_;

        autoPtr<volScalarField::Internal> wTrans_;

    // Functions

        // Initialisation

        //- Set cloud sub-models
        void setModels();

        void buildConstProps();

        void buildCellOccupancy();

        void evolveCloud(solidParcel::trackingData& td);

        //- Reset state of cloud
        void cloudReset(solidParcelCloud& c);
        
        //- Disallow default bitwise copy construct
        solidParcelCloud(const solidParcelCloud&);
        
        //- Disallow default bitwise assignment
        void operator=(const solidParcelCloud&);

public:
    // Constructors

        //- Construct given name and mesh
        solidParcelCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& Ma,
            const volScalarField& gamma,
            const volScalarField& omega,
            const dimensionedVector& g,
            const rho2ReactionThermo& thermo,
            bool readFields = true
        );
        
        //- construct given name, mesh and initialisation dictionary.
        solidParcelCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& Ma,
            const volScalarField& gamma,
            const volScalarField& omega,
            const dimensionedVector& g,
            const rho2ReactionThermo& thermo,
            const IOdictionary& solidInitialiseDict,
            const bool& clearFields
        );
    
    //- Destructor
    virtual ~solidParcelCloud();


    // Member Functions

        // Access
            labelList getTypeIDs(const dictionary& dict) const;

            //- Switch to specify if particles of the cloud can return
            //  non-zero wall distance values - true for kinematic parcels
            bool hasWallImpactDistance() const;

        // References to the mesh and databases

            //- Return the cloud type
            inline const word& cloudName() const;

            inline const label& interiorInteractionType() const;

            //- Return refernce to the mesh
            inline const fvMesh& mesh() const;

        // References to the specific data

            //- Return particle properties dictionary
            inline const IOdictionary& particleProperties() const;

            //- Return const access to the solution properties
            inline const cloudSolution& solution() const;

            //- Return reference to the sub-models dictionary
            inline const dictionary& subModelProperties() const;
            
            inline const List<word>& typeIdList() const;

            inline List<DynamicList<solidParcel*> >& cellOccupancy();

                //- Return the number of real particles represented by one parcel
                inline scalar nParticle() const;

                //- Return the radius of the domain for axisymmetric cases
                inline scalar radialExtent() const;
                
                //- Return the RWF for axisymmetric cases
                inline scalar maxRWF() const;

                //- Return the bool to identify if this is an axisymmetric simulation
                //- Measurements, etc dealt with differently if true
                inline bool axisymmetric() const;

                //- Return the number of real particles represented by one parcel
                //inline const volScalarField& nParticles() const;

                //- Return the number of real particles represented by one parcel for
                //  a given cell
                //inline scalar nParticles(const label celli) const;

                //- Return the number of real particles represented by one parcel for
                //  a given patch/face couple
                //inline scalar nParticles(const label patchi, const label facei) const;

            //- Return all of the constant properties
            inline const List<solidParcel::constantProperties>& constProps() const;
            // inline const solidParcel::constantProperties& constProps() const;

            //- Return the constant properties of the given typeId
            inline const solidParcel::constantProperties& constProps(label typeId) const;

            //- Return const access to thermo package
            inline const rho2ReactionThermo& thermo() const;

            inline const boundaryMeasurements& boundaryFluxMeasurements() const;
            inline boundaryMeasurements& boundaryFluxMeasurements(); 

            inline const solidFieldProperties& fieldProps() const;
            inline solidFieldProperties& fieldProps();

            inline const diffusionBasedSmoothing& DBS() const;
            inline diffusionBasedSmoothing& DBS();

            inline const diffusionBasedGasFiltering& DBGF() const;
            inline diffusionBasedGasFiltering& DBGF();

            inline const PackingModel& packingModel() const;
            inline PackingModel& packingModel();

            inline const DampingModel& dampingModel() const;
            inline DampingModel& dampingModel();

            inline const IsotropyModel& isotropyModel() const;
            inline IsotropyModel& isotropyModel();

            inline const HeatTransferModel& heatTransfer() const;
            inline HeatTransferModel& heatTransfer();

            //- Return reference to the random object
            inline Random& rndGen();

            //- Return the cell length scale
            inline const scalarField& cellLengthScale() const;

            inline volScalarField& theta();

            // References to the gas fields
            inline const volVectorField& UFilter() const;
            inline const volScalarField& rhoFilter() const;
            inline const volScalarField& muFilter() const;
            inline const volScalarField& TtrFilter() const;
            inline const volScalarField& CpFilter() const;
            inline const volScalarField& kappaFilter() const;
            inline const volScalarField& MaFilter() const;
            inline const volScalarField& gammaFilter() const;
            inline const volScalarField& omegaFilter() const;

        // Environmental properties

            //- Gravity
            inline const dimensionedVector& g() const;

        //- Sub-models
            inline const forceList& forces() const;
            inline forceList& forces();

            inline const solidBoundaries& boundaries() const;
            inline solidBoundaries& boundaries();

            //- Return reference to coordinate system model
            //inline const coordinateSystemType& coordSystem() const;
            //inline coordinateSystemType& coordSystem();

        // Integration schemes

            //-Return reference to velocity integration
            inline const integrationScheme& UIntegrator() const;

            //-Return reference to velocity integration
            inline const integrationScheme& TIntegrator() const;

        // Modelling options

            inline bool PreduceD() const;

        //- Source

            // Momentum

                //- Return reference to momentum source
                inline volVectorField::Internal& momentumTrans();
                // inline volVectorField& momentumTrans();

                //- Return const reference to momentum source
                inline const volVectorField::Internal& momentumTrans() const;
                // inline const volVectorField& momentumTrans() const;

                //- Return coefficient for carrier phase U equation
                inline volScalarField::Internal& momentumCoeff();

                //- Return const coefficient for carrier phase U equation
                inline const volScalarField::Internal& momentumCoeff() const;

                //- Return tmp momentum source term
                inline tmp<fvVectorMatrix> SU(volVectorField& U) const;

                //- Sensible enthalpy transfer [J/kg]
                inline volScalarField::Internal& qTrans();
                // inline volScalarField& qTrans();

                //- Sensible enthalpy transfer [J/kg]
                inline const volScalarField::Internal& qTrans() const;
                // inline const volScalarField& qTrans() const;

                //- Return coefficient for carrier phase hs equation
                inline volScalarField::Internal& qCoeff();

                //- Return const coefficient for carrier phase hs equation
                inline const volScalarField::Internal& qCoeff() const;

                //- Return internal energy source term [J/kg/m3/s]
                inline tmp<fvScalarMatrix> SE(volScalarField& e) const;

                inline volScalarField::Internal& wTrans();
                inline const volScalarField::Internal& wTrans() const;
            
        
        //- check

            //- Maximum temperature
            inline scalar Tmax() const;

            //- Minimum temperature
            inline scalar Tmin() const;

            //- Total number of parcels
            inline label nParcelsInSystem() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            // Fields

                //- Return the particle number field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> np() const;

                //- update volume fraction
                inline void updateTheta();

        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                solidParcel& parcel,
                const scalar lagrangianDt
            );

            //- Reset the cloud source terms
            void resetSourceTerms();
            void smoothSourceTerms();

            void updateGasProperties();

            void updateCellOccupancy();

            void filteringGasProperties();

            //- Cloud evolution functions
            void evolve();
            
            //- Particle motion
            // void collideMotion(solidParcel::trackingData& td);
            void MPPICMotion(solidParcel::trackingData& td);

            void axisymmetricWeighting();

            void addNewParcel
            (
                const polyMesh& owner,
                const solidParcel::constantProperties& constProps,
                const vector& position,
                const vector& U0,
                const scalar RWF,
                const label celli,
                const label tetFacei,
                const label tetPti,
                const label typeId,
                const label newParcel
            );

            void addNewParcel
            (
                const polyMesh& owner,
                const solidParcel::constantProperties& constProps,
                const vector& position,
                const vector& U0,
                const scalar T,
                const scalar RWF,
                const label celli,
                const label tetFacei,
                const label tetPti,
                const label typeId,
                const label newParcel
            );
        

        //- IO

            void info();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParcelCloudI.H"
#include "solidParcelTrackingDataI.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
