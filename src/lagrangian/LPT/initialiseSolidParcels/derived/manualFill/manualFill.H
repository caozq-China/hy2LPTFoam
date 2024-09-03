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



\*---------------------------------------------------------------------------*/

#ifndef manualFill_H
#define manualFill_H

#include "allConfigurations.H"
#include "solidParcelCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcmanualFill Declaration
\*---------------------------------------------------------------------------*/

class manualFill
:
	public configuration
{

    //- Name of file containing positions data
    const word positionsFile_;

    //- Field of parcel positions
    vectorIOField positions_;

    //- List of cell labels corresponding to injector positions
    labelList injectorCells_;

    //- List of tetFace labels corresponding to injector positions
    labelList injectorTetFaces_;

    //- List of tetPt labels corresponding to injector positions
    labelList injectorTetPts_;

public:

    //- Runtime type information
    TypeName("manualFill");

    // Constructors

        //- Construct from components
    manualFill
    (
        solidParcelCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    virtual ~manualFill() = default;


    // Member Functions

    virtual void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
