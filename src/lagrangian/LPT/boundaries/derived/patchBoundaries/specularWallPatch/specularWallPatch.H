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
    specularWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature.

    Use in the boundariesDict file.

Usage

    Example of the boundary condition specification:
    \verbatim
    patchBoundaryProperties
    {
        patchName                           wall;
    }

    boundaryModel   specularWallPatch;

    specularWallPatchProperties
    {
        temperature     300;
        velocity        (0 0 0);
    }
    \endverbatim


SourceFiles
    specularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef specularWallPatch_H
#define specularWallPatch_H

#include "solidPatchBoundary.H"
#include "solidParcelCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class specularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class specularWallPatch
:
    public solidPatchBoundary
{
    // Private data

        dictionary propsDict_;

        void setProperties();


public:

    //- Runtime type information
    TypeName("specularWallPatch");

    //- Construct from components
    specularWallPatch
    (
        const polyMesh& mesh,
        solidParcelCloud& cloud,
        const dictionary& dict
    );

    // Destructor
    virtual ~specularWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control
        virtual void controlParticle
        (
            solidParcel& p,
            solidParcel::trackingData& td
        );

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
