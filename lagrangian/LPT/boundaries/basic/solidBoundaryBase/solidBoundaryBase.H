/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef solidBoundaryBase_H
#define solidBoundaryBase_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "vectorField.H"
#include "Time.H"
#include "Random.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class solidParcelCloud;

/*---------------------------------------------------------------------------*\
                      Class Declaration
\*---------------------------------------------------------------------------*/

class solidBoundaryBase
{
protected:

    // Protected data

        //- Reference to the mesh database
        //  Note: using fvMesh so that
        const fvMesh& mesh_;

        //- Reference to the cloud
        solidParcelCloud& spc_;

        //- Reference to the time database
        const Time& time_;

        //- Patch type
        const word patchType_;

        //- Properties dictionary
        dictionary boundaryDict_;

        //- Patch name
        word patchName_;

        //- Patch ID
        label patchId_;

	// Writing out properties
	bool writeInTimeDir_;
	bool writeInCase_;


public:

    //- Runtime type information
    TypeName("solidBoundaryBase");

    //- Construct from components
    solidBoundaryBase
    (
        const polyMesh& mesh,
        solidParcelCloud& cloud,
        const dictionary& dict,
        const word& patchType
    );

    //- Destructor
    virtual ~solidBoundaryBase() = default;


    // Member Functions

        virtual void calculateProperties() = 0;

        //- Create an initial configuration
        virtual void initialConfiguration() = 0;

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);

        //- Scalar field
        virtual void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const scalarField& yData
        );

        //- Vector field
        virtual void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData
        );


        // Access

            //- Return the name of the control patch
            virtual const word& patchName() const;

            //- Return the ID of the control patch
            virtual label patchId() const;


        // Output

            virtual bool writeInTimeDir() const;
            virtual bool writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
