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

Description

\*---------------------------------------------------------------------------*/

#include "solidGeneralBoundary.H"
#include "fvMesh.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "solidParcelCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidGeneralBoundary, 0);

defineRunTimeSelectionTable(solidGeneralBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidGeneralBoundary::solidGeneralBoundary
(
    const polyMesh& mesh,
    solidParcelCloud& spc,
    const dictionary& dict
)
:
    solidBoundaryBase(mesh, spc, dict, "general"),
    faces_(),
    patchSurfaceArea_(0.0),
    cells_(),
    accumulatedParcelsToInsert_()
{
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    for(label i = 0; i < patch.size(); ++i)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidGeneralBoundary> solidGeneralBoundary::New
(
    const polyMesh& mesh,
    solidParcelCloud& spc,
    const dictionary& dict
)
{
    word modelType
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting solidGeneralBoundaryModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
//         FatalIOErrorInLookup
//         (
//             dict,
//             "boundaryModel",
//             modelType,
//             *dictionaryConstructorTablePtr_
//         ) << exit(FatalIOError);
        FatalError
            << "solidGeneralBoundary::New(const dictionary&) : " << endl
            << "    unknown solidGeneralBoundary type "
            << modelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<solidGeneralBoundary>
    (
        cstrIter()(mesh, spc, dict)
    );
}

void solidGeneralBoundary::updateTime()
{}


const labelList& solidGeneralBoundary::controlPatch() const
{
    return faces_;
}

const labelList& solidGeneralBoundary::controlZone() const
{
    return cells_;
}

const vector& solidGeneralBoundary::velocity() const
{
    return velocity_;
}


vector& solidGeneralBoundary::velocity()
{
    return velocity_;
}

const scalar& solidGeneralBoundary::temperature() const
{
    return temperature_;
}


scalar& solidGeneralBoundary::temperature()
{
    return temperature_;
}


} // End namespace Foam

// ************************************************************************* //
