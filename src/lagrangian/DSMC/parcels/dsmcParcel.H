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
    Foam::dsmcParcel

Description
    DSMC parcel class

SourceFiles
    dsmcParcelI.H
    dsmcParcel.C
    dsmcParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcParcel_H
#define dsmcParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dsmcCloud;

// Forward declaration of friend functions
/* 
Keep the two lines below. 
They are used for the forward declaration of with "Ostream"
Without these two lines, there will be a Warning during compilation.
*/
class dsmcParcel;
Ostream& operator<<(Ostream&, const dsmcParcel&);//- keep this

/*---------------------------------------------------------------------------*\
                         Class dsmcParcel Declaration
\*---------------------------------------------------------------------------*/

class dsmcParcel
:
      public particle
{
public:

    //- Class to hold DSMC particle constant properties
    class constantProperties
    {
        // Private data

            //- Particle mass[kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter[m] (constant)
            scalar d_;

            //- Rotational degrees of freedom
            label rotationalDoF_;

            //- Vibrational degrees of freedom
            label vibrationalDoF_;

            //- Viscosity index
            scalar omega_;

            //- VSS scattering parameter
            scalar alpha_;

            //- Characteristic vibrational temperature
            scalarList thetaV_;

            //- Characteristic dissociation temperature
            scalarList thetaD_;

            //- Reference vibrational relaxation number
            scalarList Zref_;

            //- Reference temperature for reference vibrational relaxation
            // number
            scalarList TrefZv_;

            //- Characteristic dissociation quantum level
            labelList charDissQuantumLevel_;

            //- Ionisation temperature
            scalar ionisationTemperature_;

            //- Number of electronic energy/degeneracy levels
            label nElectronicLevels_;

            // Charge state (0 for neutrals, 1 for ions, -1 for electrons)
            label charge_;

            //- Scalar list of electronic energy levels[J]
            scalarList electronicEnergyList_;

            //- Integer list of degeneracies at each electronic level
            labelList degeneracyList_;


    public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member functions

            //- Return const access to the particle density
            inline scalar mass() const;

            //- Return const access to the particle diameter
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the rotationalDegreesOfFreedom
            inline label rotationalDoF() const;

            //- Return the vibrationallDegreesOfFreedom
            inline label vibrationalDoF() const;

            //- Return the viscosity index
            inline scalar omega() const;

            //- Return the VSS scattering parameter
            inline scalar alpha() const;

            //- Return the characteristic vibrational temperature
            inline const scalarList& thetaV() const;

            //- Return the characteristic dissociation temperature
            inline const scalarList& thetaD() const;

            //- Return the characteristic vibrational relaxation number
            inline const scalarList& Zref() const;

            //- Return the reference temperature for the characteristic
            // vibrational relaxation number
            inline const scalarList& TrefZv() const;

            //- Return the characteristic dissociation level
            inline const labelList& charDissQuantumLevel() const;

            //- Return the ionisation temperature
            inline scalar ionisationTemperature() const;

            //- Return the number of electronic energy/degeneracy levels
            inline label nElectronicLevels() const;

            //- Return the particle's charge state
            inline label charge() const;

            //- Return the number of electronic energy/degeneracy levels
            inline const scalarList& electronicEnergyList() const;

            //- Return the number of electronic energy/degeneracy levels
            inline const labelList& degeneracyList() const;
    };


    //- Class used to pass kinematic tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Private data

    public:

        // Constructors

            //- Construct from components
            trackingData(dsmcCloud& cloud)
            :
                particle::trackingData(cloud)
            {}


    };


private:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel[m/s]
            vector U_;

            //- Radial weight of parcel (for axisymmetric)
            scalar RWF_;

            //- Rotational energy of the Parcel[J]
            scalar ERot_;

            //- Electronic energy level of the Parcel
            label ELevel_;

            //- Parcel type id
            label typeId_;

            //- newly inserted parcel
            label newParcel_;

             //- Vibrational energy of the Parcel[J]
            labelList vibLevel_;


public:

    //- Runtime type information
//     TypeName("dsmcParcel");

    friend class Cloud<dsmcParcel>;


    // Constructors

        //- Construct from components
        inline dsmcParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const vector& U,
            const scalar RWF,
            const scalar ERot,
            const label ELevel,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const label newParcel,
            const labelList& vibLevel
        );

        inline dsmcParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& U,
            const scalar RWF,
            const scalar ERot,
            const label ELevel,
            const label cellI,
            const label typeId,
            const label newParcel,
            const labelList& vibLevel
        );

        //- Construct from Istream

        dsmcParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone

        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new dsmcParcel(*this));
        }


        //- Factory class to read - construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}


            autoPtr<dsmcParcel> operator()(Istream& is) const
            {
                return autoPtr<dsmcParcel>(new dsmcParcel(mesh_, is, true));
            }
        };


    // Member Functions

        // Access

            //- Return type id
            inline label typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to radial weight
            inline scalar RWF() const;

            //- Return const access to rotational energy
            inline scalar ERot() const;

            //- Return const access to electronic energy
            inline label ELevel() const;

            //- Return const access to new parcel
            inline label newParcel() const;

            //- Return const access to vibrational energy
            inline const labelList& vibLevel() const;


        // Edit

            //- Return non - const access to typeId
            inline label& typeId();

            //- Return access to velocity
            inline vector& U();

            //- Return non-const access to radial weight
            inline scalar& RWF();

            //- Return access to rotational energy
            inline scalar& ERot();

            //- Return access to electronic energy
            inline label& ELevel();

            //- Return access to new parcel
            inline label& newParcel();

            //- Return access to vibrational energy
            inline labelList& vibLevel();


        // Main calculation loop

            // Tracking

                //- Move the parcel
                bool move
                (
                    dsmcCloud& cloud,
                    trackingData& td,
                    const scalar trackTime
                );


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch - hitting functions
            // Patch interactions

            bool hitPatch
            (
                dsmcCloud& cloud,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            // processorPatch
            void hitProcessorPatch
            (
                dsmcCloud& cloud,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                dsmcCloud& cloud,
                trackingData& td
            );


            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties(const vector& separation);


        // I - O

            static void readFields(Cloud<dsmcParcel>& c);

            static void writeFields(const Cloud<dsmcParcel>& c);


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const dsmcParcel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
