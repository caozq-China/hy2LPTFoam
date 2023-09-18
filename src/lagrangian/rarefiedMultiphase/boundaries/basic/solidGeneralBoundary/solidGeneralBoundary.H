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
    dsmcSolidGeneralBoundary

Description

SourceFiles
    dsmcSolidGeneralBoundaryI.H
    dsmcSolidGeneralBoundary.C
    dsmcSolidGeneralBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidGeneralBoundary_H
#define solidGeneralBoundary_H

// #include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "solidBoundaryBase.H"
// #include "vector.H"
// #include "volFields.H"
// #include "polyMesh.H"
// #include "solidParticleCoupling.H"
// #include "timeFluxData.H"
// #include "timeDataMeas.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class solidParticleCouplingCloud;
/*---------------------------------------------------------------------------*\
                           Class dsmcSolidGeneralBoundary Declaration
\*---------------------------------------------------------------------------*/

class solidGeneralBoundary
:
    public solidBoundaryBase
{

protected:

    // Protected data

        //- faces on the boundary
        labelList faces_;
        
        //- Type IDs (species)
        labelList typeIDs_;

        // patch surface area
        scalar patchSurfaceArea_;

        //- cells corresponding to the boundary faces
        labelList cells_;
        
//         labelList phaseStates_;
        
//         scalarList CzRatios_;

        //- set all the properties below from model if required

// 	    scalar density_;
	    vector velocity_;
	    scalar temperature_;
//	    vector USolid_;
// 	    vector FSolid_;
//	    scalar TSolid_;

	    //- targeted fields 
// 	    scalarField densities_;
// 	    vectorField velocities_;
// 	    scalarField temperatures_;

	    List<scalarField> accumulatedParcelsToInsert_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("solidGeneralBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidGeneralBoundary,
            dictionary,
            (
                const polyMesh& mesh,
                solidParticleCouplingCloud& spc,
                const dictionary& dict
            ),
            (mesh, spc, dict)
        );
   
    // Constructors

        //- Construct from components
        solidGeneralBoundary
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<solidGeneralBoundary> New
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );


    // Destructor

        virtual ~solidGeneralBoundary() = default;


    // Member Functions

        virtual void updateTime();


        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control 
 		virtual void controlParcelsBeforeMove() = 0;
        virtual void controlParcelsBeforeCollisions() = 0;
        virtual void controlParcelsAfterCollisions() = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- calculate any required properties
        virtual void calculateProperties() = 0;
        
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const vector& velocity,
            const scalarList& numDen
        );
        
        //- Insert parcels
        virtual void insertParcels
        (
            const scalar& temperature,
            const vector& velocity,
//             const scalarList& phaseStates,
            const scalarList& CzRatios
        );
        
        //- Insert parcels
        virtual void insertParcels
        (
            const scalar& temperature,
            vector& velocity,
            const scalar& nozzleRadius,
            const scalar& maxOffAxisAngle,
//             const scalarList& phaseStates,
            const scalarList& CzRatios
        );


        // Access

            virtual const labelList& controlPatch() const;

            //- return the adjacent cells
            virtual const labelList& controlZone() const;
            
            virtual const vector& velocity() const;
            virtual vector& velocity();
            
            virtual const scalar& temperature() const;
            virtual scalar& temperature();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
