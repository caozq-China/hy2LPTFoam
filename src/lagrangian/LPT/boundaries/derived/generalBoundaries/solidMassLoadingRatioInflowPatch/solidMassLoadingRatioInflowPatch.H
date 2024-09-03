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

Description


SourceFiles
    solidMassLoadingRatioInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef solidMassLoadingRatioInflowPatch_H
#define solidMassLoadingRatioInflowPatch_H

#include "solidGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class solidMassLoadingRatioInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class solidMassLoadingRatioInflowPatch
:
	public solidGeneralBoundary
{

private:

    // Private data
        dictionary solidPropsDict_;

        const scalarList massLoadingRatio_;

        scalarList parcelMass_;

        List<scalarField> numberDensities_;
        
        void setProperties();

public:

    //- Runtime type information
    TypeName("solidMassLoadingRatioInflowPatch");

    // Constructors

    //- Construct from components
    solidMassLoadingRatioInflowPatch
    (
        const polyMesh& mesh,
        solidParcelCloud& spc,
        const dictionary& dict
    );

// Destructor

    virtual ~solidMassLoadingRatioInflowPatch() = default;


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
