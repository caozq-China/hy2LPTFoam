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
    solidBoundaryMeasurements

Description

SourceFiles
    solidBoundaryMeasurementsI.H
    solidBoundaryMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef solidBoundaryMeasurements_H
#define solidBoundaryMeasurements_H

// #include "vector.H"
#include "fvMesh.H"
// #include "solidParticleCoupling.H"
// #include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class solidParticleCouplingCloud;
class solidParticleCoupling;
/*---------------------------------------------------------------------------*\
                           Class solidBoundaryMeasurements Declaration
\*---------------------------------------------------------------------------*/
// class solidParticleCouplingCloud;


class solidBoundaryMeasurements
{
    // Private data

    const fvMesh& mesh_;

    solidParticleCouplingCloud& cloud_;
    
    labelList typeIdSolids_;

    // Outer list is typeId, middle is the patch, inner is the local patch face
    List< List < vectorField > > momentumBF_;
    List< List < vectorField > > UMeanBF_;
    List< List < scalarField > > rhoNBF_;
    List< List < scalarField > > rhoMBF_;
    List< List < scalarField > > linearKEBF_;
    List< List < scalarField > > volumeFractionBF_;
    List< List < scalarField > > numberFluxBF_;
    List< List < scalarField > > massFluxBF_;
    List< List < scalarField > > conductiveHeatFluxBF_;
    
    // Private Member Functions



public:

    friend class solidParticleCoupling;

    // Constructors

        //- Construct from mesh and cloud
        solidBoundaryMeasurements
        (
            const polyMesh&,
            solidParticleCouplingCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (dsmcFoam)
        solidBoundaryMeasurements
        (
            const polyMesh&,
            solidParticleCouplingCloud& cloud,
            const bool& dummy
        );

    // Destructor

        ~solidBoundaryMeasurements() = default;

    // Member Functions

        void updateFields
        (
            solidParticleCoupling& p
        );

        void clean();

        // Access
        inline const List< List < vectorField > >& momentumBF() const;
        inline const List< List < vectorField > >& UMeanBF() const;
        inline const List< List < scalarField > >& rhoNBF() const;
        inline const List< List < scalarField > >& rhoMBF() const;
        inline const List< List < scalarField > >& linearKEBF() const;
        inline const List< List < scalarField > >& volumeFractionBF() const;
        inline const List< List < scalarField > >& numberFluxBF() const;
        inline const List< List < scalarField > >& massFluxBF() const;
        inline const List< List < scalarField > >& conductiveHeatFluxBF() const;

        // Check

        // Edit
        inline List< List < vectorField > >& momentumBF();
        inline List< List < vectorField > >& UMeanBF();
        inline List< List < scalarField > >& rhoNBF();
        inline List< List < scalarField > >& rhoMBF();
        inline List< List < scalarField > >& linearKEBF();
        inline List< List < scalarField > >& volumeFractionBF();
        inline List< List < scalarField > >& numberFluxBF();
        inline List< List < scalarField > >& massFluxBF();
        inline List< List < scalarField > >& conductiveHeatFluxBF();


        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidBoundaryMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
