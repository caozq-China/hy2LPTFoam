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
    Foam::BoikoDrag

Group

Description
    Drag coefficient used in "Shock wave interaction with a cloud of particles" 
    by Boiko et al.

SourceFiles
    BoikoDrag.C

\*---------------------------------------------------------------------------*/

#ifndef BoikoDrag_H
#define BoikoDrag_H

#include "Force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class BoikoDrag Declaration
\*---------------------------------------------------------------------------*/

class BoikoDrag
:
    public Force
{
private:
    // Private Data

    // Private Member Functions
            
        //- Drag coefficient
        scalar Cd
        (
            const scalar Re,
            const scalar Ma12
        ) const;


public:

    //- Runtime type information
    TypeName("BoikoDrag");


    // Constructors

        //- Construct from mesh
        BoikoDrag
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        BoikoDrag(const BoikoDrag& df);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new BoikoDrag(*this)
            );
        }


    //- Destructor
    virtual ~BoikoDrag();


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
