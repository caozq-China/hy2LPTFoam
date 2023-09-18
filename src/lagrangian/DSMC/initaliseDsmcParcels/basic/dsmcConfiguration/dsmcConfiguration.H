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
    Foam::dsmcConfiguration

Description


SourceFiles
    dsmcConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcConfiguration_H
#define dsmcConfiguration_H

#include "Random.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

class dsmcCloud;

/*---------------------------------------------------------------------------*\
                      Class dsmcConfiguration Declaration
\*---------------------------------------------------------------------------*/

class dsmcConfiguration
{
protected:

    // Protected data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        const dictionary& dsmcInitialiseDict_;

        Random& rndGen_;

        label nParcelsAdded_;


public:

    //- Runtime type information
    TypeName("dsmcConfiguration");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcConfiguration,
        dictionary,
        (
            dsmcCloud& cloud,
            const dictionary& dict
        ),
        (cloud, dict)
    );

    //- Construct from components
    dsmcConfiguration(dsmcCloud& cloud, const dictionary& dict);


    //- Selectors
    static autoPtr<dsmcConfiguration> New
    (
        dsmcCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~dsmcConfiguration() = default;


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        virtual label nParcelsAdded() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
