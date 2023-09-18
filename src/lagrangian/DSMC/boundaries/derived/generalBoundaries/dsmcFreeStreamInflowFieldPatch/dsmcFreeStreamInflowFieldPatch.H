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
    Foam::dsmcFreeStreamInflowFieldPatch

Description

SourceFiles
    dsmcFreeStreamInflowFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcFreeStreamInflowFieldPatch_H
#define dsmcFreeStreamInflowFieldPatch_H

#include "dsmcFreeStreamInflowPatch.H"
#include "polyMesh.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class dsmcFreeStreamInflowFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcFreeStreamInflowFieldPatch
:
    public dsmcGeneralBoundary
{
    // Private data

        dictionary propsDict_;

        // Translational temperature
        scalarField inletTransT_;

        // Rotational temperature
        scalarField inletRotT_;

        // Vibrational temperature
        scalarField inletVibT_;

        // Electronic temperature
        scalarField inletElecT_;

        vectorField inletVelocities_;

        // List of (one value) mass densities corresponding to the type Ids
        List<scalarField> numberDensities_;

        // List of species
        List<scalarField> accumulatedParcelsToInsert_;

        //- Boundary temperatures
        volScalarField boundaryTransT_;

        volScalarField boundaryRotT_;

        volScalarField boundaryVibT_;

        volScalarField boundaryElecT_;

        //- Boundary velocity
        volVectorField boundaryU_;

        //- List of inlet densities (one entry for each species)
        List<autoPtr<volScalarField>> boundaryNumberDensity_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("dsmcFreeStreamInflowFieldPatch");

    //- Construct from components
    dsmcFreeStreamInflowFieldPatch
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor
    virtual ~dsmcFreeStreamInflowFieldPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParcelsBeforeMove();
        virtual void controlParcelsBeforeCollisions();
        virtual void controlParcelsAfterCollisions();

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
