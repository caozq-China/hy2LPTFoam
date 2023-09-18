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


SourceFiles
    solidCellMeasurementsI.H
    solidCellMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef solidCellMeasurements_H
#define solidCellMeasurements_H

#include "fvMesh.H"
#include "solidParticleCoupling.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class solidCellMeasurements Declaration
\*---------------------------------------------------------------------------*/

class solidCellMeasurements
{
    // Private data

        const fvMesh& mesh_;

        solidParticleCouplingCloud& spc_;

//         scalarField collisionSeparation_;

        scalarField nColls_;
        
        scalarField nCollsTotal_;


public:

    friend class solidParticleCoupling;

    // Constructors

        //- Construct from mesh and cloud
        solidCellMeasurements
        (
            const polyMesh&,
            solidParticleCouplingCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (rarefiedMultiphaseFoam)
        solidCellMeasurements
        (
            const polyMesh&,
            solidParticleCouplingCloud& cloud,
            const bool dummy
        );


    //- Destructor
    ~solidCellMeasurements() = default;

    // Member Functions

        void updateFields(solidParticleCoupling& p);

        void clean();

//         inline const scalarField& collisionSeparation() const;

//         inline scalarField& collisionSeparation();

        inline const scalarField& nColls() const;

        inline scalarField& nColls();
        
        inline const scalarField& nCollsTotal() const;

        inline scalarField& nCollsTotal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidCellMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
