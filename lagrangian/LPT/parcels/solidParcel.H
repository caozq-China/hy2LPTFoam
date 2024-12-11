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
    Foam::solidParcel

Description

SourceFiles
    solidParcelI.H
    solidParcel.C
    solidParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParcel_H
#define solidParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"
#include "interpolation.H"
#include "AveragingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParcelCloud;

class solidParcel;

Ostream& operator<<(Ostream&, const solidParcel&);//- keep this

template<class Type>
class AveragingMethod;

/*---------------------------------------------------------------------------*\
                           Class solidParcel Declaration
\*---------------------------------------------------------------------------*/

class solidParcel
:
    public particle
{
public:
    // Private member data
    class constantProperties
    {
    //- private data
        
        //- Particle density [kg/m3] (constant)
            scalar rho0_;

            scalar d0_;

        //- thermal properties

            //- Particle initial temperature [K]
            scalar T0_;

            scalar TMin_;

            scalar TMax_;

            scalar K0_;

            //- Particle specific heat capacity [J/(kg.K)]
            scalar Cp0_;    
        
    public:
        
        // constructors
        
        inline constantProperties();
        
        //constuctor from dictionary
        inline constantProperties(const dictionary& dict);

        //- Disallow default bitwise assignment
        void operator=(const constantProperties& cps)
        {
            rho0_ = cps.rho0_;
            d0_ = cps.d0_;
            T0_ = cps.T0_;
            TMin_ = cps.TMin_;
            TMax_ = cps.TMax_;
            K0_ = cps.K0_;
            Cp0_ = cps.Cp0_;
        };

        

        // Member functions

            //- Return const access to the particle density
            inline scalar rho0() const;

            inline scalar d0() const;

            //- Return const access to the particle initial temperature [K]
            inline scalar T0() const;

            //- Return const access to minimum temperature [K]
            inline scalar TMin() const;

            //- Return const access to maximum temperature [K]
            inline scalar TMax() const;

            inline scalar K0() const;

            //- Return const access to the particle specific heat capacity
            //  [J/(kg.K)]
            inline scalar Cp0() const;

    };
    
    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<solidParcelCloud>
    {
        public:

            enum trackPart
            {
                LinearTrack,
                DampingNoTrack,
                PackingNoTrack,
                CorrectTrack
            };
            
        private:

            // Interpolators for continuous phase fields

                //- Velocity interpolator
                autoPtr<interpolation<vector>> UInterp_;

                //- Density interpolator
                autoPtr<interpolation<scalar>> rhoInterp_;

                //- Dynamic viscosity interpolator
                autoPtr<interpolation<scalar>> muInterp_;

                //- Temperature field interpolator
                autoPtr<interpolation<scalar>> TInterp_;

                //- Specific heat capacity field interpolator
                autoPtr<interpolation<scalar>> CpInterp_;

                //- Thermal conductivity field interpolator
                autoPtr<interpolation<scalar>> kappaInterp_;
                
                autoPtr<interpolation<scalar>> MaInterp_;

                //- specific heat ratio interpolator
                autoPtr<interpolation<scalar>> gammaInterp_;

                //- specific viscosity index interpolator
                autoPtr<interpolation<scalar>> omegaInterp_;

            // MPPIC Averages
                //- Velocity average
                autoPtr<AveragingMethod<vector>> uAverage_;

                //- Magnitude velocity sqyuared average
                autoPtr<AveragingMethod<scalar>> uSqrAverage_;

                //- Volume average
                autoPtr<AveragingMethod<scalar>> volumeAverage_;

                //- Radius average [ volume^(1/3) ]
                autoPtr<AveragingMethod<scalar>> radiusAverage_;

                //- Density average
                autoPtr<AveragingMethod<scalar>> rhoAverage_;

                //- Frequency average
                autoPtr<AveragingMethod<scalar>> frequencyAverage_;

                //- Mass average
                autoPtr<AveragingMethod<scalar>> massAverage_;

            //- Local gravitational or other body-force acceleration
            const vector& g_;

            //- Label specifying the current part of the tracking process
            trackPart part_;

        public:

        // Constructors

            //- Construct from components
            inline trackingData
            (
                solidParcelCloud& cloud,
                trackPart part
            );

            //- Update the MPPIC averages
            inline void updateAverages(solidParcelCloud& cloud);

            // Member functions

            //- Return conat access to the interpolator for continuous
            //  phase velocity field
            inline const interpolation<vector>& UInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase density field
            inline const interpolation<scalar>& rhoInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase dynamic viscosity field
            inline const interpolation<scalar>& muInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase temperature field
            inline const interpolation<scalar>& TInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase specific heat capacity field
            inline const interpolation<scalar>& CpInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase thermal conductivity field
            inline const interpolation<scalar>& kappaInterp() const;
            
            inline const interpolation<scalar>& MaInterp() const;

            inline const interpolation<scalar>& gammaInterp() const;

            inline const interpolation<scalar>& omegaInterp() const;

            // Return const access to the gravitational acceleration vector
            inline const vector& g() const;

            //- Const access to the tracking part label
            inline trackPart part() const;

            //- Non const access to the tracking part label
            inline trackPart& part();

    };

private:
    //- protected member data

        //- Velocity
        vector U_;

        //- Velocity correction due to collisions [m/s] used in MPPIC
        vector UCorrect_;

        scalar RWF_;

        //- Diameter [m]
        scalar d_;

        //- Density [kg/m3]
        scalar rho_;

        //- Temperature [K]
        scalar T_;

        //- conductivity
        scalar K_;

        //- Specific heat capacity [J/(kg.K)]
        scalar Cp_;

        //- solid particle type id
        label typeId_;

        label newParcel_;

        label active_;

        
        
    // gas quantities

            //- Velocity [m/s]
            vector Uc_;

            //- Density [kg/m3]
            scalar rhoc_;

            //- Viscosity [Pa.s]
            scalar muc_;

            //- Temperature [K]
            scalar Tc_;

            scalar kappac_;

            //- Specific heat capacity [J/(kg.K)]
            scalar Cpc_;

            //- Ma number
            scalar Mac_;
            
            //- specific heat ratio
            scalar gammac_;

            //- gas cell viscosity index
            scalar omegac_;

        const vector calcVelocity
        (
            trackingData& td,
            const scalar dt,           // timestep
            const label celli,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar mu,           // local carrier viscosity
            const scalar mass,         // mass
            const vector& Su,          // explicit particle momentum source
            vector& dUTrans,           // momentum transfer to carrier
            scalar& dworksTrans,
            scalar& Spu                // linearised drag coefficient
        );

        scalar calcHeatTransfer
        (
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number - surface
            const scalar kappa,        // Thermal conductivity - surface
            const scalar NCpW,         // Sum of N*Cp*W of emission species
            const scalar Sh,           // explicit particle enthalpy source
            const label celli,         // owner cell
            scalar& dhsTrans,          // sensible enthalpy transfer to carrier
            scalar& Sph                // linearised heat transfer coefficient
        );
    
public:


    // Constructors    

        //- Construct from components
        inline solidParcel
        (
            const polyMesh& mesh,
            const constantProperties& constProps,
            const vector& position,
            const vector& U0,
            const scalar RWF,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const label newParcel
        );

        inline solidParcel
        (
            const polyMesh& mesh,
            const constantProperties& constProps,
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

        //- Construct from Istream
        solidParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new solidParcel(*this));
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

            autoPtr<solidParcel> operator()(Istream& is) const
            {
                return autoPtr<solidParcel>
                (
                    new solidParcel(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return velocity
            inline const vector& U() const;

            inline const vector& UCorrect() const;

            inline scalar RWF() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to density
            inline scalar rho() const;

            //- Return const access to temperature
            inline scalar T() const;

            inline scalar K() const;

            //- Return const access to specific heat capacity
            inline scalar Cp() const;

            //- Return const access to type id
            inline label typeId() const;

            inline label newParcel() const;
            
            inline label active() const;


            //- Return const access to carrier velocity [m/s]
            inline const vector& Uc() const;

            //- Return const access to carrier density [kg/m3]
            inline scalar rhoc() const;

            //- Return const access to carrier viscosity [Pa.s]
            inline scalar muc() const;

            //- Return const access to carrier temperature
            inline scalar Tc() const;

            inline scalar kappac() const;

            //- Return const access to carrier specific heat capacity
            inline scalar Cpc() const;

            inline scalar Mac() const;

            inline scalar gammac() const;

            inline scalar omegac() const;
            
        // Edit

            //- Return access to velocity
            inline vector& U();

            //- Return access to correction velocity
            inline vector& UCorrect();

            inline scalar& RWF();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to density
            inline scalar& rho();

            //- Return access to temperature
            inline scalar& T();

            inline scalar& K();

            //- Return access to specific heat capacity
            inline scalar& Cp();

            //- Return access to type id
            inline label& typeId();

            inline label& newParcel();

            inline label& active();


        // Helper functions

            //- Cell owner mass
            inline scalar massCell(const label celli) const;

            //- Particle mass
            inline scalar mass() const;

            //- Particle moment of inertia around diameter axis
            inline scalar momentOfInertia() const;

            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline static scalar volume(const scalar d);

            //- Particle surface area
            inline scalar areaS() const;

            //- Surface area for given diameter
            inline static scalar areaS(const scalar d);

            //- Reynolds number
            inline scalar Re
            (
                const vector& U,        // particle velocity
                const scalar d,         // particle diameter
                const scalar rhoc,      // carrier density
                const scalar muc        // carrier dynamic viscosity
            ) const;

            //- Weber number
            inline scalar We
            (
                const vector& U,        // particle velocity
                const scalar d,         // particle diameter
                const scalar rhoc,      // carrier density
                const scalar sigma      // particle surface tension
            ) const;

            //- Eotvos number
            inline scalar Eo
            (
                const vector& a,        // acceleration
                const scalar d,         // particle diameter
                const scalar sigma      // particle surface tension
            ) const;

        // Main calculation loop

            //- Set cell values
            void setCellValues
            (
                trackingData& td,
                const scalar dt,
                const label celli
            );

            //- Correct cell values using latest transfer information
            void cellValueSourceCorrection
            (
                trackingData& td,
                const scalar dt,
                const label celli
            );

            //- Calculate surface thermo properties
            void calcSurfaceValues
            (
                trackingData& td,
                const scalar T,
                scalar& Ts,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappas,
                const label celli
            ) const;

            //- Update parcel properties over the time interval
            void calc
            (
                trackingData& td,
                const scalar dt,
                const label celli
            );

        // Tracking
            // scalar trackToFace(const vector& endPosition, trackingData& td);
            //- Move the parcel
            bool move(trackingData& td, const scalar trackTime);

            void LinearMove(solidParcel& p, trackingData& td, const scalar trackTime);

        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label,
                const scalar,
                const tetIndices&
            );

            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchI
            );
            
            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch& wpp,
                trackingData& td,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch& pp,
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
            virtual scalar wallImpactDistance(const vector& n) const;
            


    // I-O

        static void readFields(Cloud<solidParcel>& c);

        static void writeFields(const Cloud<solidParcel>& c);

// };
    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidParcel&);


};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
