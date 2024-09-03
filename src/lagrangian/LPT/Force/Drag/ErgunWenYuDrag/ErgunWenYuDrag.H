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
    Foam::ErgunWenYuDrag

Group

Description
    Drag coefficient used in "Shock wave interaction with a cloud of particles" 
    by Boiko et al.

SourceFiles
    ErgunWenYuDrag.C

\*---------------------------------------------------------------------------*/

#ifndef ErgunWenYuDrag_H
#define ErgunWenYuDrag_H

#include "Force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class ErgunWenYuDrag Declaration
\*---------------------------------------------------------------------------*/

class ErgunWenYuDrag
:
    public Force
{
private:
    // Private Data
        const volScalarField& thetaC_;

    // Private Member Functions
            
        //- Drag coefficient
        scalar CdRe
        (
            const scalar Re
        ) const;


public:

    //- Runtime type information
    TypeName("ErgunWenYuDrag");


    // Constructors

        //- Construct from mesh
        ErgunWenYuDrag
        (
            solidParcelCloud& cloud,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        ErgunWenYuDrag(const ErgunWenYuDrag& df);

        //- Construct and return a clone
        virtual autoPtr<Force> clone() const
        {
            return autoPtr<Force>
            (
                new ErgunWenYuDrag(*this)
            );
        }


    //- Destructor
    virtual ~ErgunWenYuDrag();


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
