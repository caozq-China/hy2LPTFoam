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
    dsmcSolidFreeStreamInflowPatch

Description

    Fully deletion boundary patch: atomisticMolecules get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    solidFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef solidFreeStreamInflowPatch_H
#define solidFreeStreamInflowPatch_H

#include "solidGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidFreeStreamInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class solidFreeStreamInflowPatch
:
	public solidGeneralBoundary
{

private:

    // Private data
        dictionary solidPropsDict_;
        
        scalarList CzRatios_;
        
//         labelList phaseStates_;
        
        void setProperties();

public:

    //- Runtime type information
    TypeName("solidFreeStreamInflowPatch");

    // Constructors

    //- Construct from components
    solidFreeStreamInflowPatch
    (
        const polyMesh& mesh,
        solidParticleCouplingCloud& spc,
        const dictionary& dict
    );

// Destructor

    virtual ~solidFreeStreamInflowPatch() = default;


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
