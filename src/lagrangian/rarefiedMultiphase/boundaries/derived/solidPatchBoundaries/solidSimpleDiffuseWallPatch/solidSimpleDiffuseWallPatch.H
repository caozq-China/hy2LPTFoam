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
    solidSimpleDiffuseWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    solidSimpleDiffuseWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef solidSimpleDiffuseWallPatch_H
#define solidSimpleDiffuseWallPatch_H

#include "solidPatchBoundary.H"
#include "solidParticleCouplingCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidSimpleDiffuseWallPatch Declaration
\*---------------------------------------------------------------------------*/

class solidSimpleDiffuseWallPatch
:
	public solidPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;

        scalar coeffRestitution_;
        
        scalar coeffFriction_;
        
        void setProperties();



public:

    //- Runtime type information
    TypeName("solidSimpleDiffuseWallPatch");

    // Constructors

        //- Construct from components
        solidSimpleDiffuseWallPatch
        (
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        virtual ~solidSimpleDiffuseWallPatch() = default;


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
