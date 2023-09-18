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

#ifndef xParticleCloud_H
#define xParticleCloud_H

#include "CloudWithModels.H"
//#include "Cloud.H"

#include "xParticle.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"

// #include "solidBoundaries.H"
// #include "solidAllConfigurations.H"
// #include "solidFieldProperties.H"

#include "dsmcCloud.H"
// #include "InterphaseCoupling.H"

// #include "solidBoundaryMeasurements.H"
// #include "solidCellMeasurements.H"
#include "Tensor.H"
#include "InteractionLists.H"
#include "CollisionModel.H"
#include "subCycleTime.H"
// #include "solidCollisionDetection.H"
// #include "solidBinaryCollisionModel.H"
// #include "solidControllers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class xParticleCloud Declaration
\*---------------------------------------------------------------------------*/
class xParticleCloud
:
      public CloudWithModels<xParticle>
{
public:

    typedef xParticleCloud baseCloudType;
    
private:
        
    // Private data
        List<word> typeIdList_;
        
        //- number of actual solid particles represented by a simulator
        scalar nRealParticles_;

        //- input for il_
        scalar maxDistance_;

        List<xParticle::constantProperties> constProps_;

        InteractionLists<xParticle> il_;

        //- Collision model
        autoPtr<CollisionModel<xParticleCloud>> collisionModel_;


        dsmcCloud* dsmcCloudPtr_;
        
        Random rndGen_;
        
        
    // Private Member Functions
        void buildConstProps();

        void moveCollide
        (
            xParticleCloud& cloud,
            xParticle::trackingData& td,
            const scalar deltaT
        );
        
        //- Disallow default bitwise copy construct
        xParticleCloud(const xParticleCloud&)= delete;
        
        //- Disallow default bitwise assignment
        void operator=(const xParticleCloud&)= delete;

public:
    
    // Constructors

        //- Construct given name and mesh, will read particles, DSMC parcels
        //- and fields from files
        xParticleCloud
        (
            const Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            dsmcCloud* dsmcCloudPtr,
            bool readFields = true
        );
        
        //- construct given name, mesh and initialisation dictionary.
        xParticleCloud
        (
            const Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& initialiseDict,
            dsmcCloud* dsmcCloudPtr,
            const bool& clearFields
        );
        
    virtual ~xParticleCloud() = default;
    // Member Functions

        // Access
            labelList getTypeIDs(const dictionary& dict) const;
            
            inline const List<word>& typeIdList() const;
            
            inline scalar nRealParticles() const;

            //- Return all of the constant properties
            inline const List<xParticle::constantProperties>& constProps() const;

            //- Return the constant properties of the given typeId
            inline const xParticle::constantProperties& constProps(label typeId) const;

            inline const InteractionLists<xParticle>& il() const;

            inline const CollisionModel<xParticleCloud>& collision() const;

            inline CollisionModel<xParticleCloud>& collision();
            
            //- Return reference to the random object 
            inline Random& rndGen();
            
        //- Sub-models
            //- Return the dsmcCloud 
            inline dsmcCloud* dsmcCloudPtr() const;

            //- clear all existing solid particles in the field
            inline void clearParticles();


            //- Cloud evolution functions
            void evolve();

            void motion
            (
                xParticleCloud& cloud,
                xParticle::trackingData& td
            );

            void calculateForce();

            void axisymmetricWeighting();
            
            void addNewParticle
            (
                const vector& position,
                const label cellI,
                const label typeId,
                const label newParticle,
                const scalar T,
                const scalar RWF,
                const vector& U,
                const vector& f,
                const vector& angularMomentum,
                const vector& torque
            );
            
            void info() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "xParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
