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
    Foam::CliftGauvinDrag

Group
    grpLagrangianIntermediateForceSubModels

Description

     Clift and Gauvin - "Motion of particles in turbulent gas streams"

SourceFiles
    CliftGauvinDrag.C

\*---------------------------------------------------------------------------*/

#ifndef CliftGauvinDrag_H
#define CliftGauvinDrag_H

#include "Force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class CliftGauvinDrag Declaration
\*---------------------------------------------------------------------------*/

class CliftGauvinDrag
:
    public Force
{
private:
    // Private Data

    // Private Member Functions
            
        //- Drag coefficient
        scalar Cd
        (
            const scalar Re
        ) const;


public:

    //- Runtime type information
    TypeName("CliftGauvinDrag");


    // Constructors

        //- Construct from mesh
        CliftGauvinDrag
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        CliftGauvinDrag(const CliftGauvinDrag& df);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new CliftGauvinDrag(*this)
            );
        }


    //- Destructor
    virtual ~CliftGauvinDrag();


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
