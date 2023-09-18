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
    Foam::dsmcBoundaries

Description
    Stores all the information for the boundary models.

    Models are defined within the boundariesDict. This class selects and builds
    the models automatically. We assume that boundaries are divided into three:

    1. Patch models:
       - actions on instantaneous collisions with a patch boundary

    2. cyclic boundary models:
       - actions on instantaneous collisions with a cyclic boundary

    3. General boundary models:
       - time-dependent actions on molecules close to the boundary


\*---------------------------------------------------------------------------*/

#ifndef dsmcBoundaries_H
#define dsmcBoundaries_H

#include "List.H"

#include "dsmcCyclicBoundary.H"
#include "dsmcPatchBoundary.H"
#include "dsmcGeneralBoundary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class dsmcBoundaries Declaration
\*---------------------------------------------------------------------------*/

class dsmcBoundaries
{
    // Private data

        const Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary dsmcBoundariesDict_;

        //- Number of boundary models
        label nPatchBoundaryModels_;
        label nCyclicBoundaryModels_;
        label nGeneralBoundaryModels_;

        //- Patch models
        PtrList<entry> patchBoundaryList_;
        List<word> patchBoundaryNames_;
        List<label> patchBoundaryIds_;
        List<fileName> pBFixedPathNames_;
        List<autoPtr<dsmcPatchBoundary>> patchBoundaryModels_;
        List<label> patchToModelId_;

        //- Cyclic boundary models
        PtrList<entry> cyclicBoundaryList_;
        List<word> cyclicBoundaryNames_;
        List<label> cyclicBoundaryIds_;
        List<fileName> cMFixedPathNames_;
        List<autoPtr<dsmcCyclicBoundary>> cyclicBoundaryModels_;
        List<label> cyclicBoundaryToModelId_;

        //- General boundary models
        PtrList<entry> generalBoundaryList_;
        List<word> generalBoundaryNames_;
        List<label> generalBoundaryIds_;
        List<fileName> gMFixedPathNames_;
        List<autoPtr<dsmcGeneralBoundary>> generalBoundaryModels_;


public:

    //- Default name of the boundaries dictionary
    static const word dictName;


    // Constructors

        //- Null Constructor
        dsmcBoundaries(const Time& t, const polyMesh& mesh);

        //- Constructor
        dsmcBoundaries(const Time& t, const polyMesh& mesh, dsmcCloud& cloud);

        ~dsmcBoundaries() = default;

        void setInitialConfig();

        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        void calculateProps();

        //- Update time information
        void updateTimeInfo();

        void outputResults();

        void checkPatchBoundaryModels(const polyMesh& mesh);
        void checkCyclicBoundaryModels(const polyMesh& mesh);


    // Member Functions

        label nPatchBoundaryModels() const;
        label nCyclicBoundaryModels() const;
        label nGeneralBoundaryModels() const;

        inline List<autoPtr<dsmcPatchBoundary>>& patchBoundaryModels();
        inline const List<autoPtr<dsmcPatchBoundary>>&
            patchBoundaryModels() const;
        inline const List<word>& patchBoundaryNames() const;
        inline const List<label>& patchBoundaryIds() const;
        inline const List<label>& patchToModelIds() const;

        inline List<autoPtr<dsmcCyclicBoundary>>& cyclicBoundaryModels();
        inline const List<autoPtr<dsmcCyclicBoundary>>&
            cyclicBoundaryModels() const;
        inline const List<word>& cyclicBoundaryNames() const;
        inline const List<label>& cyclicBoundaryIds() const;
        inline const List<label>& cyclicBoundaryToModelIds() const;


        inline List<autoPtr<dsmcGeneralBoundary>>& generalBoundaryModels();
        inline const List<autoPtr<dsmcGeneralBoundary>>&
            generalBoundaryModels() const;
        inline const List<word>& generalBoundaryNames() const;
        inline const List<label>& generalBoundaryIds() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcBoundariesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
