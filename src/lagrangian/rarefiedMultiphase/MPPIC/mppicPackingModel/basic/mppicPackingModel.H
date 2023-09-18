/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    An abstract class for selecting collisions in a DSMC simulation

Description

SourceFiles
    collisionPartnerSelectionI.H
    collisionPartnerSelection.C
    collisionPartnerSelectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef mppicPackingModel_H
#define mppicPackingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "solidParticleCouplingCloud.H"
#include "mppicParticleStressModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template <class Type>
class mppicAveragingMethod;
// class mppicParticleStressModel;
/*---------------------------------------------------------------------------*\
                           Class mppicPackingModel Declaration
\*---------------------------------------------------------------------------*/

class mppicPackingModel
{

protected:
// private:

//         const fvMesh& mesh_;
//         const dictionary& dict_;
//         word regionName_;
        
//         label regionID_;

        solidParticleCouplingCloud& spc_;
        
        //- Particle stress model
        autoPtr<mppicParticleStressModel> mppicParticleStressModel_;
        
//         autoPtr<mppicCorrectionLimitingMethods> mppicCorrectionLimitingMethods_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("mppicPackingModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            mppicPackingModel,
            dictionary,
            (
//                 const polyMesh& mesh,
                solidParticleCouplingCloud& spc,
                const dictionary& dict
            ),
            (spc, dict)
        );

    // Static data members

    // Constructors

        mppicPackingModel
        (
//             const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );
        
        //- Construct copy
//         mppicPackingModel(const mppicPackingModel& cm);

        //- Construct and return a clone
//         virtual autoPtr<mppicPackingModel> clone() const = 0;

    // Selectors

        static autoPtr<mppicPackingModel> New
        (
//             const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );

    // Destructor

        virtual ~mppicPackingModel();

    // Member Functions

//         virtual void cloudVelocityUpdateAndCorrection() = 0;
        virtual vector velocityUpdateAndCorrection
        (
            solidParticleCoupling& p,
            const scalar deltaT
        ) const = 0;
        
        virtual bool active() const =0;

        virtual void cacheFields(const bool store);
        
        virtual void calUcorrect();
        
//         virtual void evolveMPPIC(solidParticleCouplingCloud& spc);

        // Access
//         virtual const word& regionName() const;
        
        //- Return the control zone cells
//         virtual const labelList& controlZone() const;
        
        inline const mppicParticleStressModel& mppicParticleStress() const;
        
//         inline const mppicCorrectionLimitingMethods& mppicCorrectionLimiting() const;
        // Check

        // Edit
        inline mppicParticleStressModel& mppicParticleStress();
        
//         inline mppicCorrectionLimitingMethods& mppicCorrectionLimiting();

        // Write

            //- output data

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "mppicPackingModel.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
