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
    Foam::dsmcPatchBoundary

Description
    Used for applying instantaneous actions on dsmcParcels colliding with a
    boundary patch.

SourceFiles
    dsmcPatchBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcPatchBoundary_H
#define dsmcPatchBoundary_H

#include "autoPtr.H"
#include "dsmcBoundaryBase.H"
#include "dsmcParcel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class dsmcCloud;
class Time;

/*---------------------------------------------------------------------------*\
                      Class dsmcPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class dsmcPatchBoundary
:
    public dsmcBoundaryBase
{
protected:

    // Protected data

        //- Faces on the boundary
        labelList faces_;
        scalar patchSurfaceArea_;
        scalar totalPatchSurfaceArea_;

        //- Cells corresponding to the boundary faces
        labelList cells_;

        //- Set all the properties below from model if required

        scalar density_;
        vector velocity_;
        scalar temperature_;

        bool measurePropertiesAtWall_;

        //- Measured fields at the patch
        scalar preIE_;
        vector preIMom_;


public:

    //- Runtime type information
    TypeName("dsmcPatchBoundary");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcPatchBoundary,
        dictionary,
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    dsmcPatchBoundary
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );


    // Selectors

        static autoPtr<dsmcPatchBoundary> New
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dsmcPatchBoundary() = default;


    // Member Functions

        virtual void measurePropertiesBeforeControl(const dsmcParcel& p);
        virtual void measurePropertiesAfterControl
        (
            const dsmcParcel& p,
            const scalar heatOfReaction
        );

        //- Perform a diffuse wall reflection
        virtual void diffuseReflection
        (
            dsmcParcel& p,
            const scalar wallTemp,
            const vector& wallVelocity
        );

        //- perform a specular reflection
        virtual void specularReflection(dsmcParcel& p);

        //- Calculate any required properties
        virtual void calculateProperties() = 0;

        virtual void initialConfiguration() = 0;

        //- Control the dsmcParcel
        virtual void controlParticle
        (
            dsmcParcel& mol,
            dsmcParcel::trackingData& td
        ) = 0;

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;


        virtual const labelList& controlPatch() const;

        //- Return the adjacent cells
        virtual const labelList& controlZone() const;


        // Return the targeted fields

            virtual scalar density() const;
            virtual scalar& density();

            virtual const vector& velocity() const;
            virtual vector& velocity();

            virtual scalar temperature() const;
            virtual scalar& temperature();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
