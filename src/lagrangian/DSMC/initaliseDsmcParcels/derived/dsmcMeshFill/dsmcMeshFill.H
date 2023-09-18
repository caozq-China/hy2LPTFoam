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
    dsmcMeshFill

Description


SourceFiles
    dsmcMeshFillI.H
    dsmcMeshFill.C
    dsmcMeshFillIO.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcMeshFill_H
#define dsmcMeshFill_H

#include "dsmcAllConfigurations.H"
#include "vector.H"
#include "polyMesh.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class dsmcMeshFill Declaration
\*---------------------------------------------------------------------------*/

class dsmcMeshFill
:
    public dsmcConfiguration
{
public:

    //- Runtime type information
    TypeName("dsmcMeshFill");

    //- Construct from components
    dsmcMeshFill(dsmcCloud& cloud, const dictionary& dict);

    //- Destructor
    virtual ~dsmcMeshFill() = default;


    // Member Functions

        virtual void setInitialConfiguration();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
