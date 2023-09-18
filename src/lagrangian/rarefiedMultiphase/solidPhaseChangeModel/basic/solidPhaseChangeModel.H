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

#ifndef solidPhaseChangeModel_H
#define solidPhaseChangeModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
// #include "Time.H"
// #include "fvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParticleCouplingCloud;
class solidParticleCoupling;
/*---------------------------------------------------------------------------*\
                           Class collisionPartnerSelection Declaration
\*---------------------------------------------------------------------------*/

class solidPhaseChangeModel
{
    protected:

//         const fvMesh& mesh_;

        solidParticleCouplingCloud& spc_;

        //Random& rndGenS_;


    private:


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("solidPhaseChangeModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidPhaseChangeModel,
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
        //- Construct null from owner
        solidPhaseChangeModel(solidParticleCouplingCloud& owner);

        solidPhaseChangeModel
        (
//             const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );

    

    // Destructor

        virtual ~solidPhaseChangeModel() = default;
        
        // Selectors

        static autoPtr<solidPhaseChangeModel> New
        (
//             const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );

    // Member Functions
            virtual bool active() const = 0;

            virtual void initialConfiguration() = 0;

            virtual void temperatureCorrection
            (
                solidParticleCoupling& pSolid,
                const scalar& QdeltaC,
                const scalar& deltaT
            ) = 0;

            virtual void phaseCheck(solidParticleCoupling& pSolid) = 0;
            
        // Access

        // Check

        // Edit

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

// #include "collisionPartnerSelectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
