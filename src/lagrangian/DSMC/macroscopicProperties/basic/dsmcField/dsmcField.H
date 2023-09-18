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
    dsmcField

Description

SourceFiles
    dsmcFieldI.H
    dsmcField.C
    dsmcFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcField_H
#define dsmcField_H

#include "autoPtr.H"
#include "vector.H"
#include "Random.H"
#include "dsmcParcel.H"
#include "timeDataMeas.H"
#include "writeTimeData.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class dsmcField Declaration
\*---------------------------------------------------------------------------*/

class dsmcField
{
protected:

    // Protected data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        //- Time properties dictionary
        dictionary timeDict_;

        timeDataMeas timeVel_;

        fileName casePath_;

        fileName timePath_;


public:

    //- Runtime type information
    TypeName("dsmcField");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcField,
        dictionary,
        (
            const Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        ),
        (t, mesh, cloud, dict)
    );


    // Constructors

        //- Construct from components
        dsmcField
        (
            const Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcField> New
        (
            const Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~dsmcField() = default;


    // Member Functions

        void updateTime();

        const fileName& casePath() const;
        fileName& casePath();

        const fileName& timePath() const;
        fileName& timePath();

        //- Create initial field
        virtual void createField() = 0;

        //- Calculate field
        virtual void calculateField() = 0;

        //- Write field
        virtual void writeField() = 0;

        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
