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
    Foam::solidParticleCoupling

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    solidParticleCouplingI.H
    solidParticleCoupling.C
    solidParticleCouplingIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticleCoupling_H
#define solidParticleCoupling_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

// #include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParticleCouplingCloud;

template<class Type>
class mppicAveragingMethod;

/* 
Keep the two lines below. 
They are used for the forward declaration of with "Ostream"
Without these two lines, there will be a Warning during compilation.
*/
class solidParticleCoupling;
Ostream& operator<<(Ostream&, const solidParticleCoupling&);//- keep this
/*---------------------------------------------------------------------------*\
                           Class solidParticleCoupling Declaration
\*---------------------------------------------------------------------------*/

class solidParticleCoupling
:
    public particle
{
public:
    // Private member data
    class constantProperties
    {
    //- private data
        
        //- solid reference diamater in constantProperties
        scalar d_;
        
        //- Specular reflection fraction relative to all reflections
        scalar epsilon_;
        
        //- particle isothermal reflection fraction relative to total diffuse reflections
        scalar alpha_;
        
        //- material density of solid particle
        scalar rho_;
        
        //- specific heat of solid particle
        scalar Cp_;
        
        //- particle viscosity
        scalar mu_;
        
        //- particle thermal accommodation coefficient
        scalar tau_;
        
        //- particle thermal conductivity, W/(m*K)
        scalar k_;
        
        //- particle poisson ratio
        scalar poissonRatio_;
        
        //- particle elastic moduli
        scalar elasticModuli_;
        
        //- particle equilibrium melting temperature
        scalar Tm_;
        
        //- particle maximum temperature for homogeneous crystallization
        scalar Tf_;
        
        //- nonspherical solid particle volume
        //- if this is 0, then the solid particle is spherical
        //- else, the solid particle is nonsperical
        scalar nonSphericalVolume_;
        
        //- nonspherical solid particle superfical area
        //- if this is 0, then the solid particle is spherical
        //- else, the solid particle is nonsperical
        scalar nonSphericalSuperficalArea_;
        
        scalar effD_;
        
        scalar effRho_;
        
        scalar massSphere_;
        
        bool nonSperical_;
        
        
        
    public:
        
        // constructors
        
        inline constantProperties();
        
        //constuctor from dictionary
        inline constantProperties(const dictionary& dict);

        // Member functions
        inline scalar d() const;
        
        inline scalar epsilon() const;
        
        inline scalar alpha() const;
        
        inline scalar rho() const;
        
        inline scalar Cp() const;
        
        inline scalar mu() const;
        
        inline scalar tau() const;
        
        inline scalar k() const;
        
        inline scalar poissonRatio() const;
        
        inline scalar elasticModuli() const;
        
        inline scalar Tm() const;
        
        inline scalar Tf() const;
        
        inline scalar sigmaT() const;
        
        inline scalar massSphere() const;
        
        inline bool nonSperical() const;
        
    };
    
    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        public:
            enum trackPart
            {
                LinearTrack,
//                 DampingNoTrack,
//                 PackingNoTrack,
                CorrectTrack
            };
            
            
            
            trackingData(solidParticleCouplingCloud& spc)
            :
                  particle::trackingData(spc)
            {}
            
            //- Const access to the tracking part label
            inline trackPart part() const;

            //- Non const access to the tracking part label
            inline trackPart& part();
            
        private:
            //- Label specifying the current part of the tracking process
            trackPart part_;
    };
    

    //- used when MPPIC model is activated
    template<class CloudType>
    class TrackingData
    :
        public particle::trackingData
    {
    public:
            // Constructors
            inline TrackingData(CloudType& cloud);
            
            //- Update the MPPIC averages
            inline void updateMPPICAverages(CloudType& cloud);

    private:
        // MPPIC Averages

            //- Volume average
            autoPtr<mppicAveragingMethod<scalar> > volAve_;
            
            //- Radius average [ volume^(1/3) ]
            autoPtr<mppicAveragingMethod<scalar> > radiusAve_;
            
            //- Density average
            autoPtr<mppicAveragingMethod<scalar> > rhoAve_;

            //- Velocity average
            autoPtr<mppicAveragingMethod<vector> > uAve_;

            //- Magnitude velocity sqyuared average
            autoPtr<mppicAveragingMethod<scalar> > uSqrAve_;
            
            //- Frequency average
            autoPtr<mppicAveragingMethod<scalar> > frequencyAve_;

            //- Mass average
            autoPtr<mppicAveragingMethod<scalar> > massAve_;

    };

private:
    //- protected member data
    
    	//- Velocity
        vector U_;

        //- angular velocity
        vector omega_;
        
        //- Velocity correction due to collisions [m/s] in MPPIC
        vector UCorrect_;
        
        //- Aerodynamic Force addition on the solid particle of each time step
        vector F_;
        
        //- instantaneous diameter of an individual solid particle, 
        //- e.g. particle phase change may cause the change in particle size
        scalar D_;
        
        //- Particle temperature
        scalar T_;
        
        //- Radial weight of parcel (for axisymmetric)
        scalar RWF_;
        
        //- ratio of the crystallization front radius to particle radius
        scalar CzRatio_;

        //- numerical sequence
        label numSeq_;
        
        //- solid particle type id
        label typeID_;
        
        //- phase state of the solid particles
        //- 0: pure solid phase 1: solidification 2: melting 3: pure liquid phase
        label phaseState_;
        
        //- adding new solid particles
        label newSolidParticle_;
        
        //- build an object for translating solid diameter to wallImpactDistance
        //- in solidParticleCoupling.C
//         constantProperties constProperties_;
        
    
public:

    friend class Cloud<solidParticleCoupling>;
        
        //- Construct from components
        inline solidParticleCoupling
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const vector& U,
            const vector& omega,
            const vector& UCorrect,
            const vector& F,
            const scalar D,
            const scalar T,
            const scalar RWF,
            const scalar CzRatio,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label numSeq,
            const label typeID,
            const label phaseState,
            const label newSolidParticle
        );
        
        inline solidParticleCoupling
        (
            const polyMesh& mesh,
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
            const label typeID,
            const label phaseState,
            const label newSolidParticle
        );

        //- Construct from Istream
        solidParticleCoupling
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new solidParticleCoupling(*this));
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

            autoPtr<solidParticleCoupling> operator()(Istream& is) const
            {
                return autoPtr<solidParticleCoupling>
                (
                    new solidParticleCoupling(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access
            //- Return velocity
//             inline const vector& UCorrect() const;
            
            inline const vector& U() const;

            inline const vector& omega() const;
            
            inline const vector& UCorrect() const;

            inline const vector& F() const;
            
            //- Return temperature
            inline scalar T() const;
            
            inline scalar D() const;
            
            inline scalar RWF() const;
            
            inline scalar CzRatio() const;

            inline label numSeq() const;
			
            inline label typeID() const;
            
            inline label phaseState() const;
            
            inline label newSolidParticle() const;
            
        // Edit

            //- Return velocity
//             inline vector& UCorrect();
            
            inline vector& U();

            inline vector& omega();
            
            inline vector& UCorrect();
			
            inline vector& F();
            
            //- Return temperature
            inline scalar& T();
            
            inline scalar& D();
            
            inline scalar& RWF();
            
            inline scalar& CzRatio();

            inline label& numSeq();
            
            inline label& typeID();
            
            inline label& phaseState();

            inline label& newSolidParticle();
            
            
        
        // Tracking

            //- Move
            bool move
            (
                solidParticleCouplingCloud& spc,
                trackingData& td, 
                const scalar trackTime
            );

        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                solidParticleCouplingCloud& spc,
                trackingData& td
            );
            
            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                solidParticleCouplingCloud& spc,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                solidParticleCouplingCloud& spc,
                trackingData& td
            );


            //- Overridable function to handle the particle hitting a polyPatch
//             void hitPatch
//             (
//                 const polyPatch&,
//                 trackingData& td
//             );

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

        static void readFields(Cloud<solidParticleCoupling>& c);

        static void writeFields(const Cloud<solidParticleCoupling>& c);

// };
    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidParticleCoupling&);


};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleCouplingI.H"

#include "solidParticleCouplingMPPICTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
