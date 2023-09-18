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

#ifndef absoluteMethod_H
#define absoluteMethod_H

#include "mppicCorrectionLimitingMethods.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class absoluteMethod Declaration
\*---------------------------------------------------------------------------*/

class absoluteMethod
:
    public mppicCorrectionLimitingMethods
{
    // Private data
    
        dictionary propsDict_;
        
        //- Coefficient of restitution
        scalar e_;


public:

    //- Runtime type information
    TypeName("absoluteMethod");


    //- Constructors

        //- Construct from components
        absoluteMethod
        (
            const dictionary& dict
        );
        
//         absoluteMethod(const absoluteMethod& cl);

        //- Construct and return a clone
        virtual autoPtr<mppicCorrectionLimitingMethods> clone() const
        {
            return autoPtr<mppicCorrectionLimitingMethods>
            (
                new absoluteMethod(*this)
            );
        }


    //- Destructor
    virtual ~absoluteMethod();


    //- Member Functions

        vector limitVelocity
        (
            const vector uP,
            const vector dU,
            const vector uMean
        ) const;
        

        //- Collision stress derivaive w.r.t. the volume fraction
//         tmp<Field<scalar> > dTaudTheta
//         (
//             const Field<scalar>& alpha,
//             const Field<scalar>& rho,
//             const Field<scalar>& uRms
//         ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
