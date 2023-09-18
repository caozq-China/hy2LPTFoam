/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::dsmcGeneralBoundary

Description

SourceFiles
    dsmcGeneralBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcGeneralBoundary_H
#define dsmcGeneralBoundary_H

#include "autoPtr.H"
#include "dsmcBoundaryBase.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class dsmcCloud;

/*---------------------------------------------------------------------------*\
                     Class dsmcGeneralBoundary Declaration
\*---------------------------------------------------------------------------*/

class dsmcGeneralBoundary
:
    public dsmcBoundaryBase
{
protected:

    // Protected data

        //- Faces on the boundary
        labelList faces_;

        //- Type IDs (species)
        labelList typeIds_;

        //- Patch surface area
        scalar patchSurfaceArea_;

        //- Cells corresponding to the boundary faces
        labelList cells_;

        // Set all the properties below from model if required
        scalar density_;
        vector velocity_;
        scalar translationalTemperature_;
        scalar rotationalTemperature_;
        scalar vibrationalTemperature_;
        scalar electronicTemperature_;
        List<scalarField> accumulatedParcelsToInsert_;


public:

    //- Runtime type information
    TypeName("dsmcGeneralBoundary");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcGeneralBoundary,
        dictionary,
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    dsmcGeneralBoundary
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );


    //- Selectors
    static autoPtr<dsmcGeneralBoundary> New
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~dsmcGeneralBoundary() = default;


    // Member Functions

        virtual void updateTime();

        //- Create an initial configuration
        virtual void initialConfiguration() = 0;

        //- Control
        virtual void controlParcelsBeforeMove() = 0;
        virtual void controlParcelsBeforeCollisions() = 0;
        virtual void controlParcelsAfterCollisions() = 0;

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- Calculate any required properties
        virtual void calculateProperties() = 0;

        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalar& transT,
            const vector& velocity,
            const scalarList& numDen
        );

        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalarField& transT,
            const vectorField& velocity,
            const List<scalarField>& numDen
        );
        
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalar& transT,
            const vectorField& velocity,
            const scalar& numDen,
            const scalarField& molFractions
        );
        
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalarField& transT,
            const vectorField& velocity,
            const scalarField& numDen,
            const scalarField& molFractions
        );
        
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalar& transT,
            const vectorField& velocity,
            const List<scalarField>& numDen
        );

        //- Insert parcels
        virtual void insertParcels
        (
            const scalar& transT,
            const scalar& rotT,
            const scalar& vibT,
            const scalar& elecT,
            const vector& velocity
        );

        //- Insert parcels
        virtual void insertParcels
        (
            const scalarField& transT,
            const scalarField& rotT,
            const scalarField& vibT,
            const scalarField& elecT,
            const vectorField& velocity
        );
        
        //- Insert parcels
        virtual void insertParcels
        (
            const scalar& transT,
            const vectorField& velocity
        );
        
        //- Insert parcels
        virtual void insertParcels
        (
            const scalarField& transT,
            const vectorField& velocity
        );

        // Access

            virtual const labelList& controlPatch() const;

            //- Return the adjacent cells
            virtual const labelList& controlZone() const;


            // Return the targeted fields

            virtual scalar density() const;
            virtual scalar& density();

            virtual const vector& velocity() const;
            virtual vector& velocity();

            virtual scalar translationalTemperature() const;
            virtual scalar& translationalTemperature();

            virtual scalar rotationalTemperature() const;
            virtual scalar& rotationalTemperature();

            virtual scalar vibrationalTemperature() const;
            virtual scalar& vibrationalTemperature();

            virtual scalar electronicTemperature() const;
            virtual scalar& electronicTemperature();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
