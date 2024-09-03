/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

SourceFiles
    solidCyclicBoundary.H
    solidCyclicBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef solidCyclicBoundary_H
#define solidCyclicBoundary_H

#include "autoPtr.H"
#include "solidBoundaryBase.H"
#include "runTimeSelectionTables.H"
#include "solidParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidCyclicBoundary Declaration
\*---------------------------------------------------------------------------*/

class solidCyclicBoundary
:
    public solidBoundaryBase
{

protected:

    // Protected data

        //- corresponding/neighbouring coupled patch 
        word neighbPatchName_;
        label neighbPatchId_;

        labelList faces_;

        labelList coupledFacesA_;
        labelList coupledFacesB_;

        labelList cellsA_;
        labelList cellsB_;

    // Private Member Functions

        void getCoupledFaces(const polyMesh& mesh);

public:

    //- Runtime type information
    TypeName("solidCyclicBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidCyclicBoundary,
            dictionary,
            (
                const polyMesh& mesh,
                solidParcelCloud& spc,
                const dictionary& dict
            ),
            (mesh, spc, dict)
        );
   
    // Constructors

        //- Construct from components
        solidCyclicBoundary
        (
            const polyMesh& mesh,
            solidParcelCloud& spc,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<solidCyclicBoundary> New
        (
            const polyMesh& mesh,
            solidParcelCloud& spc,
            const dictionary& dict
        );


    // Destructor

        virtual ~solidCyclicBoundary() = default;


    // Member Functions

        virtual void calculateProperties() = 0;

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control the molecule
        virtual void controlParticle
        (
            solidParcel& p,
            solidParcel::trackingData& td
        ) = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        // Access

            //- retrun faces on patch
            virtual const labelList& controlPatch() const;


            //- return the adjacent cells
            virtual const labelList& controlZone() const;

            virtual const labelList& allFaces() const;
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
