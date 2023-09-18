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
    dsmcSolidSpecularWallPatch

Description

     Specular wall boundary patch: reflecting molecules back in by reversing their
        normal velocity vector to the boundary

SourceFiles
    dsmcSolidSpecularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef solidSpecularWallPatch_H
#define solidSpecularWallPatch_H

#include "solidPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
// #include "volFields.H"
#include "solidParticleCouplingCloud.H"
// #include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidSpecularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class solidSpecularWallPatch
:
	public virtual solidPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;


public:

    //- Runtime type information
    TypeName("solidSpecularWallPatch");

    // Constructors

    //- Construct from components
    solidSpecularWallPatch
	(
        const polyMesh& mesh,
        solidParticleCouplingCloud& cloud,
        const dictionary& dict
	);

    // Destructor

 		virtual ~solidSpecularWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle(solidParticleCoupling& p, solidParticleCoupling::trackingData& td);

        virtual void injectParticlesFromWall();

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
