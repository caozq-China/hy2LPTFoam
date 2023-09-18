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
    solidField

Description

SourceFiles
    solidFieldI.H
    solidField.C
    solidFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidField_H
#define solidField_H

// #include "IOdictionary.H"
// #include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
// #include "volFields.H"

#include "Random.H"
#include "solidParticleCoupling.H"
#include "timeDataMeas.H"
// #include "propertyField.H"
#include "writeTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidField Declaration
\*---------------------------------------------------------------------------*/

class solidField
{

protected:

    // Protected data

        const fvMesh& mesh_;

        solidParticleCouplingCloud& cloud_;

        //- subDictionary containing the time properties
        dictionary timeDict_;

        timeDataMeas timeVel_;

        fileName casePath_;

        fileName timePath_;


    // Private Member Functions

public:

    //- Runtime type information
    TypeName("solidField");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidField,
            dictionary,
            (
                const Time& t,
                const polyMesh& mesh,
                solidParticleCouplingCloud& cloud,
                const dictionary& dict
            ),
            (t, mesh, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        solidField
        (
            const Time& t,
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<solidField> New
        (
            const Time& t,
            const polyMesh& mesh,
            solidParticleCouplingCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        virtual ~solidField() = default;


    // Member Functions


        void updateTime();

        const fileName& casePath() const;
        fileName& casePath();

        const fileName& timePath() const;
        fileName& timePath();

        //- create initial field
        virtual void createField() = 0;

        //- calculate field
        virtual void calculateField() = 0;

        //- write field
        virtual void writeField() = 0;

//         void updateBasicFieldProperties(const dictionary&);

        virtual void updateProperties(const dictionary& dict);
        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
