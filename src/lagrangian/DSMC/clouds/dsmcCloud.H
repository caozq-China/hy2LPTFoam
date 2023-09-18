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

Class
    Foam::dsmcCloud

Description
    Base class for DSMC cloud.

SourceFiles
    dsmcCloudI.H
    dsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCloud_H
#define dsmcCloud_H

#include "CloudWithModels.H"

#include "dsmcParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "BinaryCollisionModel.H"

#include "dsmcBoundaries.H"
#include "dsmcAllConfigurations.H"
#include "dsmcFieldProperties.H"
#include "dsmcFaceTracker.H"
#include "dsmcControllers.H"

#include "dsmcReactions.H"
#include "collisionPartnerSelection.H"
#include "boundaryMeasurements.H"
#include "cellMeasurements.H"

#include "CloudFunctionObjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class dsmcCloud Declaration
\*---------------------------------------------------------------------------*/


class dsmcCloud
:
    public CloudWithModels<dsmcParcel>
{
public:

    typedef dsmcCloud baseCloudType;

private:

    // Private data

        //- A list of unique instances of molecule types in the
        //- simulation.
        //  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (N2 O2 CO2) N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Number of real atoms/molecules represented by a parcel
        scalar nParticle_;


        // Axisymmetric set-ups

            //- Bool to identify is this is an axisymmetric simulation
            //  Movements dealt with differently if true
            bool axisymmetric_;

            //- Radial extent of the axisymmetric domain
            scalar radialExtent_;

            //- Maximum radial scaling factor for axisymmetric simulations
            scalar maxRWF_;


        //- Write out the information to terminal every nTerminalOutputs
        label nTerminalOutputs_;

        //- A data structure holding electron number density
        scalarList rhoNMeanElectron_;

        //- A data structure holding electron mass density
        scalarList rhoMMeanElectron_;

        //- A data structure holding all species mass density
        scalarList rhoMMean_;

        //- A data structure holding electron momentum
        List<vector> momentumMeanElectron_;

        //- A data structure holding all species momentum
        List<vector> momentumMean_;

        //- A data structure holding electron linearKE
        scalarList linearKEMeanElectron_;

        //- A data structure holding electron temperature
        scalarList electronTemperature_;

        //- A data structure holding cell velocity
        List<vector> cellVelocity_;

        //- A field holding the value of (sigmaT * cR)max for each
        //- cell (see Bird p220).
        //  Initialised with the parcels, updated as required, and read in on
        //  start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Parcel constant properties - one for each type
        List<dsmcParcel::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //- Controllers
        dsmcControllers controllers_;

        //- Fields
        dsmcFieldProperties fields_;

        //- Boundaries
        dsmcBoundaries boundaries_;

        //- Fluxes - face tracking
        dsmcFaceTracker trackingInfo_;


        // References to the cloud sub - models

            //- Binary collision model
            autoPtr<BinaryCollisionModel> binaryCollisionModel_;

            //- Collision selection routine model
            autoPtr<collisionPartnerSelection> collisionPartnerSelectionModel_;

            //- Reactions
            dsmcReactions reactions_;

            //- Boundary measurements
            boundaryMeasurements boundaryMeas_;

            //- Used to calculate mean collision separation
            cellMeasurements cellMeas_;

            CloudFunctionObjectList<baseCloudType> functions_;


    // Private Member Functions


        label pickFromCandidateList(DynamicList<label>& candidatesInCell);

        // whenever you pick from the candidate list
        // always update the candidate sublist
        void updateCandidateSubList
        (
            const label candidate,
            DynamicList<label>& candidatesInSubCell
        );

        label pickFromCandidateSubList
        (
            DynamicList<label>& candidatesInCell,
            DynamicList<label>& candidatesInSubCell
        );

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Remove electrons after adding their velocities to a DynamicList
        void removeElectrons();

        //- Add an electron where an ion is found, after the move function
        void addElectrons();

        //- Calculate collisions between molecules
        void collisions();

        //- Disallow default bitwise copy construct
        dsmcCloud(const dsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const dsmcCloud&);
        
        //- for rarefiedMultiphase
        friend class solidParticleCouplingCloud;

public:

    // Constructors

        //- Construct given name and mesh, will read Parcels and fields from
        //- file
        dsmcCloud
        (
            const Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary
        dsmcCloud
        (
            const Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& dsmcInitialiseDict,
            const bool clearFields
        );

    //- Destructor
    virtual ~dsmcCloud() = default;


    // Member Functions

        //- Return the particle IDs corresponding to the particle names read
        //- from the input dictionary
        labelList getTypeIDs(const dictionary& dict) const;

        //- Return the idList
        inline const List<word>& typeIdList() const;

        //- Return the number of real particles represented by one parcel
        inline scalar nParticle() const;

        //- Return the radius of the domain for axisymmetric cases
        inline scalar radialExtent() const;

        //- Return the radial weighting factor for axisymmetric cases
        inline scalar maxRWF() const;

        //- Return the radial weighting factor for axisymmetric cases
        //  Note: assumes radius given by the y compoent!
        inline scalar axiRWF(const point& p) const;

        //- Return the bool to identify if this is an axisymmetric simulation
        //  Measurements, etc dealt with differently if true
        inline bool axisymmetric() const;

        //- Return the number of timesteps to write out to terminal at
        inline label nTerminalOutputs() const;

        //- Return the sigmaTcRMax field.  non - const access to allow
        //- updating.
        inline volScalarField& sigmaTcRMax();

        //- Return the collision selection remainder field.  non - const
        //- access to allow updating.
        inline scalarField& collisionSelectionRemainder();

        //- Return all of the constant properties
        inline const List<dsmcParcel::constantProperties>& constProps() const;

        //- Return the constant properties of the given typeId
        inline const dsmcParcel::constantProperties& constProps
        (
            const label typeId
        ) const;

        //- Return refernce to the random object
        inline Random& rndGen();

        // Kinetic theory helper functions

        //- Generate a random velocity sampled from the Maxwellian speed
        //- distribution
        vector equipartitionLinearVelocity
        (
            const scalar temperature,
            const scalar mass
        );

        //- Generate a random velocity sampled from the Chapman-Enskog speed
        // distribution
        vector chapmanEnskogLinearVelocity
        (
            scalar maxMacroscopic,
            scalar pressure,
            scalar temperature,
            scalar mass,
            vector heatflux,
            tensor stress
        );

        //- Generate a random rotational energy, sampled from the
        //- equilibrium distribution (Bird eqn 11.22 and 11.23 and
        //- adapting code from DSMC3.FOR)
        scalar equipartitionRotationalEnergy
        (
            const scalar temperature,
            const scalar rotationalDoF
        );

        //- Generate a random vibrational energy
        labelList equipartitionVibrationalEnergyLevel
        (
            const scalar temperature,
            const scalar vibrationalDoF,
            const label typeId
        );

        //- Generate a random electronic energy level
        label equipartitionElectronicLevel
        (
            const scalar temperature,
            const List<label>& degeneracyList,
            const List<scalar>& electronicEnergyList,
            const label typeId
        );

        //- Generate a LB post - collision rotational energy
        scalar postCollisionRotationalEnergy
        (
            const scalar rotationalDof,
            const scalar ChiB
        );

        //- Generate a LB post - collision vibrational energy
        label postCollisionVibrationalEnergyLevel
        (
            const bool postReaction,
            const label vibLevel,
            const label iMax,
            const scalar thetaV,
            const scalar thetaD,
            const scalar refTempZv,
            const scalar omega,
            const scalar Zref,
            const scalar Ec
        );

        //- Generate a LB post-collision electronic energy level
        label postCollisionElectronicEnergyLevel
        (
            const scalar Ec,
            const label jMax,
            const scalar omega,
            const List<scalar>& EElist,
            const List<label>& gList
        );

        //- Average particle speed from the Maxwellian distribution
        inline scalar maxwellianAverageSpeed
        (
            const scalar temperature,
            const scalar mass
        ) const;

        inline tmp<scalarField> maxwellianAverageSpeed
        (
            const scalarField& temperature,
            const scalar mass
        ) const;

        //- RMS particle speed
        inline scalar maxwellianRMSSpeed
        (
            const scalar temperature,
            const scalar mass
        ) const;

        inline tmp<scalarField> maxwellianRMSSpeed
        (
            const scalarField& temperature,
            const scalar mass
        ) const;

        //- Most probable speed
        inline scalar maxwellianMostProbableSpeed
        (
            const scalar temperature,
            const scalar mass
        ) const;

        inline tmp<scalarField> maxwellianMostProbableSpeed
        (
            const scalarField& temperature,
            const scalar mass
        ) const;


    // Sub-models

        //- Return reference to binary elastic collision model
        inline const BinaryCollisionModel& binaryCollision() const;

        //- Return non-const reference to binary elastic collision model
        inline BinaryCollisionModel& binaryCollision();

        inline const dsmcBoundaries& boundaries() const;
        inline dsmcBoundaries& boundaries();

        inline const dsmcFieldProperties& fields() const;
        inline dsmcFieldProperties& fields();

        inline const dsmcControllers& controllers() const;
        inline dsmcControllers& controllers();

        inline const dsmcFaceTracker& tracker() const;
        inline dsmcFaceTracker& tracker();

        inline const dsmcReactions& reactions() const;
        inline dsmcReactions& reactions();

        inline const boundaryMeasurements& boundaryFluxMeasurements() const;
        inline boundaryMeasurements& boundaryFluxMeasurements();

        inline const cellMeasurements& cellPropMeasurements() const;
        inline cellMeasurements& cellPropMeasurements();


        // Check

            //- Total mass in system
            inline scalarList infoMeasurements() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;

            void axisymmetricWeighting();


        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const scalar RWF,
                const scalar ERot,
                const label ELevel,
                const label cellI,
                const label typeId,
                const label newParcel,
                const labelList vibLevel
            );

            //- Redistribute vibrational energy after a reaction
            scalar energyRatio
            (
                const scalar ChiA,
                const scalar ChiB
            );

            //- Redistribute vibrational energy after a reaction
            scalar PSIm
            (
                const scalar DOFm,
                const scalar DOFtot
            );

            //- Evolve the cloud (move, collide)
            void evolve();

            //- Retrun how often to write to terminal
            // (used for the solver to access this)
            label nTerminalOutputs();

            //- Clear the Cloud
            inline void clear();


            CloudFunctionObjectList<baseCloudType>& functions()
            {
                return functions_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
