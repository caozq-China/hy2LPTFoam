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
    boundaryMeasurements

Description

SourceFiles
    boundaryMeasurementsI.H
    boundaryMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryMeasurements_H
#define boundaryMeasurements_H

#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class solidParcelCloud;
class solidParcel;
/*---------------------------------------------------------------------------*\
                           Class boundaryMeasurements Declaration
\*---------------------------------------------------------------------------*/


class boundaryMeasurements
{
    // Private data

    const fvMesh& mesh_;

    solidParcelCloud& cloud_;
    
    labelList typeIds_;

    // Outer list is typeId, middle is the patch, inner is the local patch face
    List< List < vectorField > > momentumBF_;
    List< List < vectorField > > UMeanBF_;
    List< List < vectorField > > fDBF_;
    List< List < scalarField > > rhoNBF_;
    List< List < scalarField > > rhoMBF_;
    List< List < scalarField > > qBF_;
    
    // Private Member Functions



public:

    friend class solidParcel;

    // Constructors

        //- Construct from mesh and cloud
        boundaryMeasurements
        (
            const polyMesh&,
            solidParcelCloud& cloud
        );

        //- Construct from mesh, cloud and boolean
        boundaryMeasurements
        (
            const polyMesh&,
            solidParcelCloud& cloud,
            const bool& dummy
        );

    // Destructor

        ~boundaryMeasurements() = default;

    // Member Functions

        void updateFields
        (
            solidParcel& p
        );

        void clean();

        // Access
        inline const List< List < vectorField > >& momentumBF() const;
        inline const List< List < vectorField > >& UMeanBF() const;
        inline const List< List < vectorField > >& fDBF() const;
        inline const List< List < scalarField > >& rhoNBF() const;
        inline const List< List < scalarField > >& rhoMBF() const;
        inline const List< List < scalarField > >& qBF() const;

        // Check

        // Edit
        inline List< List < vectorField > >& momentumBF();
        inline List< List < vectorField > >& UMeanBF();
        inline List< List < vectorField > >& fDBF();
        inline List< List < scalarField > >& rhoNBF();
        inline List< List < scalarField > >& rhoMBF();
        inline List< List < scalarField > >& qBF();


        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
