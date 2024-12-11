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


\*---------------------------------------------------------------------------*/

#ifndef solidPatchBoundary_H
#define solidPatchBoundary_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "solidParcel.H"
#include "solidBoundaryBase.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
class polyMesh;
class fvMesh;
class solidParcelCloud;
class Time;

/*---------------------------------------------------------------------------*\
                           Class Declaration
\*---------------------------------------------------------------------------*/

class solidPatchBoundary
:
    public solidBoundaryBase
{
protected:

    // Protected data
        
        labelList faces_;
        scalar patchSurfaceArea_;
        scalar totalPatchSurfaceArea_;
        
        //- cells corresponding to the boundary faces
        labelList cells_;
        
        //- set all the properties below from model if required
        vector velocity_;
        scalar density_;
        scalar temperature_;
        
        bool measurePropertiesAtWall_;
        bool calculateHeatConduction_;
        
        scalar preIkineticE_;
        vector preIMom_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("solidPatchBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidPatchBoundary,
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
        solidPatchBoundary
        (
            const polyMesh& mesh,
            solidParcelCloud& spc,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<solidPatchBoundary> New
        (
            const polyMesh& mesh,
            solidParcelCloud& spc,
            const dictionary& dict
        );


    // Destructor

        virtual ~solidPatchBoundary() = default;


    // Member Functions

        virtual void measurePropertiesBeforeControl(const solidParcel& p);
        virtual void measurePropertiesAfterControl(const solidParcel& p);

        //- calculate any required properties
        virtual void calculateProperties() = 0;

        virtual void initialConfiguration() = 0;

        //- control the solid particle
        virtual void controlParticle
        (
            solidParcel& p,
            solidParcel::trackingData& td
        ) = 0;

//         virtual void injectParticlesFromWall() = 0;
	
        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        // Access

        virtual const labelList& controlPatch() const;

        //- return the adjacent cells
        virtual const labelList& controlZone() const;

        //- return the targeted fields 
        virtual const scalar& density() const;
        virtual scalar& density();

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
