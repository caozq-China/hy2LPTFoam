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
    Foam::HendersonDrag

Group

Description

SourceFiles
    HendersonDrag.C

\*---------------------------------------------------------------------------*/

#ifndef HendersonDrag_H
#define HendersonDrag_H

#include "Force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class HendersonDrag Declaration
\*---------------------------------------------------------------------------*/

class HendersonDrag
:
    public Force
{
private:
    // Private Data

    // Private Member Functions
        scalar molecularSpeedRatio(const scalar Map, const scalar gamma) const;

        scalar TpOverTc(const scalar Tp, const scalar Tc) const;

        scalar volumeFractionCorrectionFactor(const scalar thetaCCell) const;
            
        //- Drag coefficient
        scalar Cd
        (
            const scalar Ma,
            const scalar Re,
            const scalar gamma,//- ratio of specific heats of gas
            const scalar Tc,   //- gas phase temperature
            const scalar Tp
        ) const;


public:

    //- Runtime type information
    TypeName("HendersonDrag");


    // Constructors

        //- Construct from mesh
        HendersonDrag
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        HendersonDrag(const HendersonDrag& df);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new HendersonDrag(*this)
            );
        }


    //- Destructor
    virtual ~HendersonDrag();


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

#endif

// ************************************************************************* //
