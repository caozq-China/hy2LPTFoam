/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::cellMeasurements

Description

    This class measures the number of collisions in each cell and the distance
    between colliding particles. This is used to return the mean collision
    separation, a key parameter for ensuring good DSMC practice.


    Reference:
        G. A. Bird. "Visual DSMC Program for Two-Dimensional and Axially
        Symmetric Flows", The DS2V Program User’s Guide, Version 3.8,
        GAB Consulting, Sydney, Australia, October 2006.

SourceFiles
    cellMeasurementsI.H
    cellMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef cellMeasurements_H
#define cellMeasurements_H

#include "fvMesh.H"
#include "dsmcParcel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cellMeasurements Declaration
\*---------------------------------------------------------------------------*/

class cellMeasurements
{
    // Private data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        scalarField collisionSeparation_;

        scalarField nColls_;


public:

    friend class dsmcParcel;

    // Constructors

        //- Construct from mesh and cloud
        cellMeasurements
        (
            const polyMesh&,
            dsmcCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (dsmcFoam)
        cellMeasurements
        (
            const polyMesh&,
            dsmcCloud& cloud,
            const bool dummy
        );


    //- Destructor
    ~cellMeasurements() = default;

    // Member Functions

        void updateFields(dsmcParcel& p);

        void clean();

        inline const scalarField& collisionSeparation() const;

        inline scalarField& collisionSeparation();

        inline const scalarField& nColls() const;

        inline scalarField& nColls();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
