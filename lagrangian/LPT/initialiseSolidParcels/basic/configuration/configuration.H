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

#ifndef configuration_H
#define configuration_H

#include "runTimeSelectionTables.H"
#include "Random.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

class solidParcelCloud;
    
    
    
/*---------------------------------------------------------------------------*\
                           Class configuration Declaration
\*---------------------------------------------------------------------------*/

class configuration
{

protected:

    // Protected data

        const fvMesh& mesh_;

        solidParcelCloud& cloud_;

        const dictionary& initialiseDict_;

        label nParticlesAdded_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("configuration");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            configuration,
            dictionary,
            (
                solidParcelCloud& cloud,
                const dictionary& dict
            ),
            (cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        configuration
        (
            solidParcelCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<configuration> New
        (
            solidParcelCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~configuration() = default;


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        const label& nParticlesAdded() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
