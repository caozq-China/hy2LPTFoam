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
    dsmcFreeStreamInflowPatch

Description

    This boundary condition only suitable for axisymmetric case.

SourceFiles
    dsmcFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcMorrisNozzleExitPatch_H
#define dsmcMorrisNozzleExitPatch_H

#include "dsmcGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcMorrisNozzleExitPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcMorrisNozzleExitPatch
:
	public dsmcGeneralBoundary
{

private:

    // Private data
        dictionary propsDict_;
        
        // Boundary temperatures
        vectorField velocities_;
        scalarField transTemp_;
        scalarField rotatTemp_;
        scalarField vibrTemp_;
        scalarField elecTemp_;
        scalar nozzleExitRadius_;

        // list of number densities corresponding to the type Ids
        List< scalarField> numberDensities_;
        
        
        scalarList rhoX_;
        scalarList rhoy_;
        
        scalarList TX_;
        scalarList TY_;
        
        scalarList axialUx_;
        scalarList axialUy_;
        
        scalarList radialVx_;
        scalarList radialVy_;
        
        
        scalar linearInterpolation(scalarList& listX,scalarList& listY, scalar Rratio);
        
        void initialiseRhoX();
        void initialiseRhoY();
        void initialiseTX();
        void initialiseTY();
        void initialiseUx();
        void initialiseUy();
        void initialiseVx();
        void initialiseVy();
        
        void setProperties();

public:

    //- Runtime type information
    TypeName("dsmcMorrisNozzleExitPatch");

    // Constructors

    //- Construct from components
    dsmcMorrisNozzleExitPatch
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

// Destructor

    virtual ~dsmcMorrisNozzleExitPatch() = default;


// Member Functions

    virtual void initialConfiguration();

    //- calculate any required properties
    virtual void calculateProperties();

    //- control the mols
    virtual void controlParcelsBeforeMove();
    virtual void controlParcelsBeforeCollisions();
    virtual void controlParcelsAfterCollisions();

    //- output data
    virtual void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    //- update properties from a modified dictionary
    virtual void updateProperties(const dictionary& dict);



        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
