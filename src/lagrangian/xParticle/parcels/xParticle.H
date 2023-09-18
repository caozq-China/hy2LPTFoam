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
    Foam::xParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    xParticleI.H
    xParticle.C
    xParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef xParticle_H
#define xParticle_H

// #include "particle.H"
#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

#include "CollisionRecordList.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"

// #include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef CollisionRecordList<vector, vector> collisionRecordList;
typedef vectorFieldCompactIOField pairDataFieldCompactIOField;
typedef vectorFieldCompactIOField wallDataFieldCompactIOField;

class xParticle;
class xParticleCloud;

// template<class Type>
// class mppicAveragingMethod;

/* 
Keep the two lines below. 
They are used for the forward declaration of with "Ostream"
Without these two lines, there will be a Warning during compilation.
*/
// class solidParticleCoupling;
Ostream& operator<<(Ostream&, const xParticle&);//- keep this
/*---------------------------------------------------------------------------*\
                           Class xParticle Declaration
\*---------------------------------------------------------------------------*/

class xParticle
:
//     public particle
    public particle
{
public:
    // Private member data
    class constantProperties
    {
    //- private data
        scalar d_;

        //- material density of solid particle
        scalar rho_;
        
        //- specific heat of solid particle
        scalar Cp_;
        
        scalar mass_;
        
        
    public:
        
        // constructors
        
        //- Null constructor, allows List of constantProperties to be
        //  created before the contents is initialised
        inline constantProperties();
        
        //constuctor from dictionary
        inline constantProperties(const dictionary& dict);

        // Member functions
        inline scalar d() const;

        inline scalar rho() const;
        
        inline scalar Cp() const;
        
        inline scalar mass() const;
        
    };
    
    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // label specifying which part of the integration algorithm is taking
        label part_;

    public:
            // Constructors
            trackingData(xParticleCloud& cloud, label part)
            :
                particle::trackingData(cloud),
                part_(part)
            {}

            //- Member functions

            inline label part() const
            {
                return part_;
            }
    };

private:
    //- protected member data
        //- solid particle type id
        label typeID_;

        //- adding new particles
        label newParticle_;
        
        //- Particle temperature
        scalar T_;
        
        //- Radial weight of parcel (for axisymmetric)
        scalar RWF_;

        vector U_;

        //- angular velocity
        vector omega_;

        //- force
        vector f_;

        //- Angular momentum of Parcel in global reference frame [kg m2/s]
        vector angularMomentum_;

        //- Torque on particle due to collisions in global
        //  reference frame [Nm]
        vector torque_;

        //- Particle collision records
        collisionRecordList collisionRecords_;
    
public:

    friend class Cloud<xParticle>;

         //- String representation of properties
        AddToPropertyList
        (
            xParticle,
            " (fx fy fz)"
          + " (angularMomentumx angularMomentumy angularMomentumz)"
          + " (torquex torquey torquez)"
          + " collisionRecordsPairAccessed"
          + " collisionRecordsPairOrigProcOfOther"
          + " collisionRecordsPairOrigIdOfOther"
          + " (collisionRecordsPairData)"
          + " collisionRecordsWallAccessed"
          + " collisionRecordsWallPRel"
          + " (collisionRecordsWallData)"
        );

        //- Construct from components
        inline xParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeID,
            const label newParticle,
            const scalar T,
            const scalar RWF,
            const vector& U,
            // const vector& omega,
            const vector& f,
            const vector& angularMomentum,
            const vector& torque
        );

        inline xParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label typeID,
            const label newParticle,
            const scalar T,
            const scalar RWF,
            const vector& U,
            // const vector& omega,
            const vector& f,
            const vector& angularMomentum,
            const vector& torque
        );

        //- Construct from Istream
        xParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new xParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<xParticle> operator()(Istream& is) const
            {
                return autoPtr<xParticle>::New(mesh_, is, true);
            }
        };


    // Member Functions

        // Access
            
            //- Return temperature
            inline label typeID() const;

            inline label newParticle() const;

            inline scalar T() const;
            
            inline scalar RWF() const;

            inline const vector& U() const;

            inline const vector& f() const;

            inline const vector& angularMomentum() const;

            inline const vector& torque() const;

            //- Return const access to the collision records
            inline const collisionRecordList& collisionRecords() const;

            inline vector omega() const;
            
        // Edit
            inline label& typeID();

            inline label& newParticle();

            //- Return temperature
            inline scalar& T();
            
            inline scalar& RWF();

            inline vector& U();

            // inline vector& omega();

            inline vector& f();

            inline vector& angularMomentum();
            
            inline vector& torque();

            //- Return access to collision records
            inline collisionRecordList& collisionRecords();
            
        
        // Tracking

            //- Move
            bool move
            (
                xParticleCloud&,
                trackingData&,
                const scalar trackTime
            );

        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                xParticleCloud& cloud,
                trackingData& td
            );
            
            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                xParticleCloud& cloud,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                xParticleCloud& cloud,
                trackingData& td
            );


            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
             virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
             virtual void transformProperties(const vector& separation);

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
//             virtual scalar wallImpactDistance(const vector& n) const;
            


    // I-O

        static void readFields(Cloud<xParticle>& c);

        static void writeFields(const Cloud<xParticle>& c);

// };
    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const xParticle&);


};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "xParticleI.H"

// #include "solidParticleCouplingMPPICTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
