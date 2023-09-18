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
    dsmcBoundaries

Description

	Stores all the information for the boundary models. Models are defined
    within the boundariesDict. This class selects & builds the models 
    automatically. We assume that boundaries are divided into three:

    1. Patch models -- actions on instantaneous collisions with a patch boundary 

    2. cyclic boundary models -- actions on instantaneous collisions with 
                                a cyclic boundary 

    3. General boundary models -- time-dependent actions on molecules close to
                                  the boundary


\*---------------------------------------------------------------------------*/

#ifndef solidBoundaries_H
#define solidBoundaries_H

#include "List.H"

#include "solidCyclicBoundary.H"
#include "solidPatchBoundary.H"
#include "solidGeneralBoundary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcBoundaries Declaration
\*---------------------------------------------------------------------------*/

class solidBoundaries
{
    // Private data


        const Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary solidBoundariesDict_;

        //- number of boundary models
        label nSolidPatchBoundaryModels_;
        label nSolidCyclicBoundaryModels_;
        label nSolidGeneralBoundaryModels_;

        //- patch models
        PtrList<entry> solidPatchBoundaryList_;
        List<word> solidPatchBoundaryNames_;
        List<label> solidPatchBoundaryIds_;
        List<fileName> solidPBFixedPathNames_;
        List< autoPtr<solidPatchBoundary> > solidPatchBoundaryModels_;
        List<label> solidPatchToModelId_;

        //- cyclic boundary models
        PtrList<entry> solidCyclicBoundaryList_;
        List<word> solidCyclicBoundaryNames_;
        List<label> solidCyclicBoundaryIds_;
        List<fileName> solidCMFixedPathNames_;
        List< autoPtr<solidCyclicBoundary> > solidCyclicBoundaryModels_;
        List<label> solidCyclicBoundaryToModelId_;

        //- solid general boundary models 
        PtrList<entry> solidGeneralBoundaryList_;
        List<word> solidGeneralBoundaryNames_;
        List<label> solidGeneralBoundaryIds_;
        List<fileName> solidGMFixedPathNames_;
        List< autoPtr<solidGeneralBoundary> > solidGeneralBoundaryModels_;





public:
    
    //- Default name of the boundaries dictionary
    static const word dictName;

    // Constructors
    
        //- Null Constructor 
        solidBoundaries
        (
            const Time& t,
            const polyMesh& mesh
        );

        //- Constructor for gnemdFOAM
        solidBoundaries
        (
            const Time& t,
            const polyMesh& mesh,
            solidParticleCouplingCloud& spc
        );

        ~solidBoundaries() = default;

        void setInitialConfig();

        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();
        
        void releaseParticlesFromWall();

        void calculateProps();

        //- update time information
        void updateTimeInfo();

        void outputResults();

        void checkSolidPatchBoundaryModels(const polyMesh& mesh);
        void checkSolidCyclicBoundaryModels(const polyMesh& mesh);


    // Member Functions

        label nSolidPatchBoundaryModels() const;
        label nSolidCyclicBoundaryModels() const;
        label nSolidGeneralBoundaryModels() const;

        inline List< autoPtr<solidPatchBoundary> >& solidPatchBoundaryModels();
        inline const List< autoPtr<solidPatchBoundary> >& solidPatchBoundaryModels() const;
        inline const List<word>& solidPatchBoundaryNames() const;
        inline const List<label>& solidPatchBoundaryIds() const;
        inline const List<label>& solidPatchToModelIds() const;

        inline List< autoPtr<solidCyclicBoundary> >& solidCyclicBoundaryModels();
        inline const List< autoPtr<solidCyclicBoundary> >& solidCyclicBoundaryModels() const;
        inline const List<word>& solidCyclicBoundaryNames() const;
        inline const List<label>& solidCyclicBoundaryIds() const;
        inline const List<label>& solidCyclicBoundaryToModelIds() const;
        
        inline List< autoPtr<solidGeneralBoundary> >& solidGeneralBoundaryModels();
        inline const List< autoPtr<solidGeneralBoundary> >& solidGeneralBoundaryModels() const;
        inline const List<word>& solidGeneralBoundaryNames() const;
        inline const List<label>& solidGeneralBoundaryIds() const;

    // Ostream Operator

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidBoundariesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
