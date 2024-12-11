/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::SinghDrag

Group

Description

SourceFiles
    SinghDrag.C

\*---------------------------------------------------------------------------*/

#ifndef SinghDrag_H
#define SinghDrag_H

#include "Force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class SinghDrag Declaration
\*---------------------------------------------------------------------------*/
class SinghDrag
:
    public Force
{
private:
    // Private Data

        //- Reference to the gas volume fraction field
        // const volScalarField& thetaC_;

        bool volumeFractionCorrection_;

    // Private Member Functions

        scalar ReInfinity(const scalar d, const vector Uc_, const scalar rhoc, const scalar muc) const;

        // scalar volumeFractionCorrectionFactor(const scalar thetaCCell) const;
            
        //- Drag coefficient
        scalar Cd
        (
            const scalar Ma,
            const scalar Re,
            const scalar gamma,//- ratio of specific heats of gas
            const scalar Tc,   //- gas phase temperature
            const scalar Tp,
            const scalar omega //- index of power law of gas viscosity
        ) const;


public:

    //- Runtime type information
    TypeName("SinghDrag");


    // Constructors

        //- Construct from mesh
        SinghDrag
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        SinghDrag(const SinghDrag& df);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new SinghDrag(*this)
            );
        }

    //- Destructor
    virtual ~SinghDrag();


    // Member Functions

        // Evaluation
            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const solidParcel& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
