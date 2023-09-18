/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::ParticleStressModels::HarrisCrighton

Description
    Inter-particle stress model of Harris and Crighton

    The stress value takes the following form:
    \f[
        \frac{P_s \alpha^\beta}{ \mathrm{max} \left( \alpha_{pack} - \alpha ,
        \epsilon ( 1 - \alpha ) \right) }
    \f]
    Here, \f$\alpha\f$ is the volume fraction of the dispersed phase, and the
    other values are modelling constants. A small value \f$\epsilon\f$ is used
    to limit the denominator to ensure numerical stability.

    Reference:
    \verbatim
        "Solitons, solitary waves, and voidage disturbances in gas-fluidized
        beds"
        S Harris and D Crighton,
        Journal of Fluid Mechanics
        Volume 266, Pages 243-276, 1994
    \endverbatim

SourceFiles
    HarrisCrighton.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitPacking_H
#define ExplicitPacking_H

#include "mppicPackingModel.H"
#include "mppicAveragingMethod.H"
#include "mppicCorrectionLimitingMethods.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class mppicAveragingMethod;

class mppicCorrectionLimitingMethods;
/*---------------------------------------------------------------------------*\
                        Class ExplicitPacking Declaration
\*---------------------------------------------------------------------------*/

class ExplicitPacking
:
    public mppicPackingModel
{
private:
    // Private data

        //- Velocity average
        const mppicAveragingMethod<vector>* velocityAve_;

        //- Volume fraction average
        const mppicAveragingMethod<scalar>* volumeAverage_;
        
        autoPtr<mppicAveragingMethod<scalar> > stressAverage_;
        
        autoPtr<mppicCorrectionLimitingMethods> mppicCorrectionLimitingMethods_;


    // Private member functions


public:

    //- Runtime type information
    TypeName("ExplicitPacking");


    //- Constructors

        //- Construct from components
        ExplicitPacking
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );


    //- Destructor
    ~ExplicitPacking();


    //- Member Functions
    
        
        
//         void cloudVelocityUpdateAndCorrection();
        //- Calculate the velocity correction
        vector velocityUpdateAndCorrection
        (
            solidParticleCoupling& p,
            const scalar deltaT
        ) const;
        
        bool active() const;
        
        //- Calculate the inter particles stresses
        void cacheFields(const bool store);
        
        void calUcorrect();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } //End namespace mppicPackingModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
