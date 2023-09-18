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

\*---------------------------------------------------------------------------*/

#include "dsmcBoundaries.H"
#include "emptyPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::dsmcBoundaries::dictName("boundariesDict");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dsmcBoundaries::dsmcBoundaries(const Time& t, const polyMesh& mesh)
:
    time_(t),
    dsmcBoundariesDict_
    (
        IOobject
        (
            dsmcBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(),
    patchBoundaryNames_(),
    patchBoundaryIds_(),
    pBFixedPathNames_(),
    patchBoundaryModels_(),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(),
    cyclicBoundaryNames_(),
    cyclicBoundaryIds_(),
    cMFixedPathNames_(),
    cyclicBoundaryModels_(),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(),
    generalBoundaryNames_(),
    generalBoundaryIds_(),
    gMFixedPathNames_(),
    generalBoundaryModels_()
{}


dsmcBoundaries::dsmcBoundaries
(
    const Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    dsmcBoundariesDict_
    (
        IOobject
        (
            dsmcBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(dsmcBoundariesDict_.lookup("dsmcPatchBoundaries")),
    patchBoundaryNames_(patchBoundaryList_.size()),
    patchBoundaryIds_(patchBoundaryList_.size()),
    pBFixedPathNames_(patchBoundaryList_.size()),
    patchBoundaryModels_(patchBoundaryList_.size()),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(dsmcBoundariesDict_.lookup("dsmcCyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(dsmcBoundariesDict_.lookup("dsmcGeneralBoundaries")),
    generalBoundaryNames_(generalBoundaryList_.size()),
    generalBoundaryIds_(generalBoundaryList_.size()),
    gMFixedPathNames_(generalBoundaryList_.size()),
    generalBoundaryModels_(generalBoundaryList_.size())
{
    Info << "Creating the boundary models: " << nl << endl;

    // Patch boundaries

    if (patchBoundaryModels_.size() > 0)
    {
        forAll(patchBoundaryModels_, p)
        {
            const entry& boundaryI = patchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();

            patchBoundaryModels_[p] = autoPtr<dsmcPatchBoundary>
            (
                dsmcPatchBoundary::New(mesh, cloud, boundaryIDict)
            );

            patchBoundaryNames_[p] = patchBoundaryModels_[p]->type();
            patchBoundaryIds_[p] = p;
            ++nPatchBoundaryModels_;
        }
    }

    checkPatchBoundaryModels(mesh);

    // Cyclic boundaries

    if (cyclicBoundaryModels_.size() > 0)
    {
        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();

            cyclicBoundaryModels_[c] = autoPtr<dsmcCyclicBoundary>
            (
                dsmcCyclicBoundary::New(mesh, cloud, boundaryIDict)
            );

            cyclicBoundaryNames_[c] = cyclicBoundaryModels_[c]->type();
            cyclicBoundaryIds_[c] = c;
            ++nCyclicBoundaryModels_;
        }
    }

    checkCyclicBoundaryModels(mesh);

    // General boundaries

    if (generalBoundaryModels_.size() > 0)
    {
        forAll(generalBoundaryModels_, g)
        {
            const entry& boundaryI = generalBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();

            generalBoundaryModels_[g] = autoPtr<dsmcGeneralBoundary>
            (
                dsmcGeneralBoundary::New(mesh, cloud, boundaryIDict)
            );

            generalBoundaryNames_[g] = generalBoundaryModels_[g]->type();
            generalBoundaryIds_[g] = g;
            ++nGeneralBoundaryModels_;
        }
    }

    // Creating directories
    if (nPatchBoundaryModels_ > 0)
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if (!isDir(boundariesPath))
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/dsmc
        fileName dsmcBoundariesPath(boundariesPath/"dsmc");

        if (!isDir(dsmcBoundariesPath))
        {
            mkDir(dsmcBoundariesPath);
        }

        // directory: case/boundaries/dsmc/patchBoundaryModels
        fileName patchBoundaryModelsPath
        (
            dsmcBoundariesPath/"patchBoundaryModels"
        );

        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);
        }

        forAll(patchBoundaryModels_, p)
        {
            if (patchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/dsmc/
                // patchBoundaryModels/<patchBoundaryModel>
                fileName patchBoundaryModelPath
                (
                    patchBoundaryModelsPath/patchBoundaryNames_[p]
                );

                if (!isDir(patchBoundaryModelPath))
                {
                    mkDir(patchBoundaryModelPath);
                }

                const word& patchName = patchBoundaryModels_[p]->patchName();

                // directory: case/controllers/dsmc/
                // patchBoundaryModels/<patchBoundaryModel>/<patchName>

                fileName patchPath(patchBoundaryModelPath/patchName);

                if (!isDir(patchPath))
                {
                    mkDir(patchPath);
                }

                pBFixedPathNames_[p] = patchPath;
            }
        }
    }


    // Creating directories
    if (nCyclicBoundaryModels_ > 0)
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if (!isDir(boundariesPath))
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/dsmc
        fileName dsmcBoundariesPath(boundariesPath/"dsmc");

        if (!isDir(dsmcBoundariesPath))
        {
            mkDir(dsmcBoundariesPath);
        }

        // directory: case/boundaries/dsmc/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath
                    (dsmcBoundariesPath/"cyclicBoundaryModels");

        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);
        }

        forAll(cyclicBoundaryModels_, c)
        {
            if (cyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/dsmc/
                // cyclicBoundaryModels/<cyclicBoundaryModel>
                fileName cyclicBoundaryModelPath
                (
                    cyclicBoundaryModelsPath/cyclicBoundaryNames_[c]
                );

                if (!isDir(cyclicBoundaryModelPath))
                {
                    mkDir(cyclicBoundaryModelPath);
                }

                const word& patchName = cyclicBoundaryModels_[c]->patchName();

                // directory: case/controllers/dsmc
                // cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>

                fileName patchPath(cyclicBoundaryModelPath/patchName);

                if (!isDir(patchPath))
                {
                    mkDir(patchPath);
                }

                cMFixedPathNames_[c] = patchPath;
            }
        }
    }

    // Creating directories
    if (nGeneralBoundaryModels_ > 0)
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if (!isDir(boundariesPath))
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/dsmc
        fileName dsmcBoundariesPath(boundariesPath/"dsmc");

        if (!isDir(dsmcBoundariesPath))
        {
            mkDir(dsmcBoundariesPath);
        }

        // directory: case/boundaries/dsmc/cyclicBoundaryModels
        fileName generalBoundaryModelsPath
        (
            dsmcBoundariesPath/"generalBoundaryModels"
        );

        if (!isDir(generalBoundaryModelsPath))
        {
            mkDir(generalBoundaryModelsPath);
        }

        forAll(generalBoundaryModels_, g)
        {
            if (generalBoundaryModels_[g]->writeInCase())
            {
                // directory: case/boundaries/dsmc
                // generalBoundaryModels/<generalBoundaryModel>
                fileName generalBoundaryModelPath
                (
                    generalBoundaryModelsPath/generalBoundaryNames_[g]
                );

                if (!isDir(generalBoundaryModelPath))
                {
                    mkDir(generalBoundaryModelPath);
                }

                const word& patchName = generalBoundaryModels_[g]->patchName();

                // directory: case/controllers/dsmc
                // generalBoundaryModels/<generalBoundaryModel>/<patchName>

                fileName patchPath(generalBoundaryModelPath/patchName);

                if (!isDir(patchPath))
                {
                    mkDir(patchPath);
                }

                gMFixedPathNames_[g] = patchPath;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcBoundaries::checkCyclicBoundaryModels(const polyMesh& mesh)
{
    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (isA<cyclicPolyPatch>(patch))
        {
            label patchIndex = patch.index();

            forAll(cyclicBoundaryModels_, c)
            {
                const label patchId = cyclicBoundaryModels_[c]->patchId();

                if (patchIndex == patchId)
                {
                    ++nPolyPatches;
                    cyclicBoundaryToModelId_[patchi] = c;
                }
            }
        }
    }

    if (nPolyPatches != nCyclicBoundaryModels_)
    {
        FatalErrorInFunction
            << " Number of cyclic boundary models = " << nCyclicBoundaryModels_
            << " chosen in the boundaryiesDict are inconsistent."
            << abort(FatalError);
    }
}


void dsmcBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)
{
    // Check that all poly-patches defined within blockMeshDict,
    // each have one model.

    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (!polyPatch::constraintType(patch.type()))
        {
            ++nPolyPatches;

            label patchIndex = patch.index();

            label nPatches = 0;

            forAll(patchBoundaryModels_, p)
            {
                const label patchId = patchBoundaryModels_[p]->patchId();

                if (patchIndex == patchId)
                {
                    ++nPatches;
                    patchToModelId_[patchi] = p;
                }
            }

            if (nPatches > 1)
            {
                FatalIOErrorInFunction(dsmcBoundariesDict_)
                    << " Only one patch boundary model per poly-patch, [name: "
                    << patch.name()
                    << "]. Number of models chosen for this patch are: "
                    << nPatches
                    << abort(FatalIOError);
            }
        }
    }

    if (nPolyPatches != nPatchBoundaryModels_)
    {
        FatalIOErrorInFunction(dsmcBoundariesDict_)
            << " Number of poly-patches = "  << nPolyPatches
            << " are not equal to the number of patch "
            << "models = " << nPatchBoundaryModels_
            << abort(FatalIOError);
    }
}


void dsmcBoundaries::updateTimeInfo()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->updateTime();
    }
}


void dsmcBoundaries::setInitialConfig()
{
    for (auto& model : patchBoundaryModels_)
    {
        model->initialConfiguration();
    }

    for (auto& model : cyclicBoundaryModels_)
    {
        model->initialConfiguration();
    }

    for (auto& model : generalBoundaryModels_)
    {
        model->initialConfiguration();
    }
}


void dsmcBoundaries::calculateProps()
{
    for (auto& model : patchBoundaryModels_)
    {
        model->calculateProperties();
    }

    for (auto& model : cyclicBoundaryModels_)
    {
        model->calculateProperties();
    }

    for (auto& model : generalBoundaryModels_)
    {
        model->calculateProperties();
    }
}


void dsmcBoundaries::controlBeforeMove()
{
    // Impose model after calculation of forces
    for (auto& model: generalBoundaryModels_)
    {
        model->controlParcelsBeforeMove();
    }
}


void dsmcBoundaries::controlBeforeCollisions()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->controlParcelsBeforeCollisions();
    }
}


void dsmcBoundaries::controlAfterCollisions()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->controlParcelsAfterCollisions();
    }
}


void dsmcBoundaries::outputResults()
{
    if (!time_.writeTime())
    {
        return;
    }

    {
        List<fileName> timePathNames(pBFixedPathNames_.size());

        if (nPatchBoundaryModels_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath
                (time_.path()/time_.timeName()/"uniform");

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            fileName boundariesTimePath(uniformTimePath/"boundaries");

            if (!isDir(boundariesTimePath))
            {
                mkDir(boundariesTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries/dsmc
            fileName dsmcTimePath(boundariesTimePath/"dsmc");

            if (!isDir(dsmcTimePath))
            {
                mkDir(dsmcTimePath);
            }

            // directory: case/<timeDir>/uniform
                    // boundaries/dsmc/patchBoundaryModels
            fileName dsmcPatchBoundaryModelsTimePath
                (dsmcTimePath/"patchBoundaryModels");

            if (!isDir(dsmcPatchBoundaryModelsTimePath))
            {
                mkDir(dsmcPatchBoundaryModelsTimePath);
            }

            forAll(patchBoundaryModels_, p)
            {
                if
                (
                    patchBoundaryModels_[p]->writeInTimeDir()
                 || patchBoundaryModels_[p]->writeInCase()
                )
                {
                    // creating directory for different zones but of the
                    // same model
                    const word& patchName =
                        patchBoundaryModels_[p]->patchName();

                    // directory: case/<timeDir>/uniform/
                                // controllers/dsmc/
                                // patchBoundaryModels/
                                //<patchBoundaryModel>/<patchName>
                    fileName patchTimePath
                    (
                        dsmcPatchBoundaryModelsTimePath/patchName
                    );

                    if (!isDir(patchTimePath))
                    {
                        mkDir(patchTimePath);
                    }

                    timePathNames[p] = patchTimePath;

                    patchBoundaryModels_[p]->output
                        (pBFixedPathNames_[p], timePathNames[p]);
                }
            }
        }
    }

    {
        List<fileName> timePathNames(cMFixedPathNames_.size());

        if (nCyclicBoundaryModels_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath
                (time_.path()/time_.timeName()/"uniform");

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            fileName boundariesTimePath(uniformTimePath/"boundaries");

            if (!isDir(boundariesTimePath))
            {
                mkDir(boundariesTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries/dsmc
            fileName dsmcTimePath(boundariesTimePath/"dsmc");

            if (!isDir(dsmcTimePath))
            {
                mkDir(dsmcTimePath);
            }

            // directory: case/<timeDir>/uniform
            // boundaries/dsmc/cyclicBoundaryModels
            fileName dsmcCyclicBoundaryModelsTimePath
                (dsmcTimePath/"cyclicBoundaryModels");

            if (!isDir(dsmcCyclicBoundaryModelsTimePath))
            {
                mkDir(dsmcCyclicBoundaryModelsTimePath);
            }

            forAll(cyclicBoundaryModels_, c)
            {
                if
                (
                    cyclicBoundaryModels_[c]->writeInTimeDir()
                 || cyclicBoundaryModels_[c]->writeInCase()
                )
                {
                    // directory: case/<timeDir>/uniform/controllers
                    // dsmc/cyclicBoundaryModels/<cyclicBoundaryModel>
                    fileName cMTimePath
                    (
                        dsmcCyclicBoundaryModelsTimePath
                       /cMFixedPathNames_[c]
                    );

                    if (!isDir(cMTimePath))
                    {
                        mkDir(cMTimePath);
                    }

                    // creating directory for different zones but
                    // of the same model
                    const word& patchName =
                        cyclicBoundaryModels_[c]->patchName();

                    // directory: case/<timeDir>/uniform/controllers
                    // dsmc/cyclicBoundaryModels
                    //<cyclicBoundaryModel>/<patchName>
                    fileName patchTimePath(cMTimePath/patchName);

                    if (!isDir(patchTimePath))
                    {
                        mkDir(patchTimePath);
                    }

                    timePathNames[c] = patchTimePath;

                    cyclicBoundaryModels_[c]->output
                    (
                        cMFixedPathNames_[c],
                        timePathNames[c]
                    );
                }
            }
        }
    }

    {
        List<fileName> timePathNames(gMFixedPathNames_.size());

        if (nGeneralBoundaryModels_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath
            (
                time_.path()/time_.timeName()/"uniform"
            );

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            fileName boundariesTimePath(uniformTimePath/"boundaries");

            if (!isDir(boundariesTimePath))
            {
                mkDir(boundariesTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries/dsmc
            fileName dsmcTimePath(boundariesTimePath/"dsmc");

            if (!isDir(dsmcTimePath))
            {
                mkDir(dsmcTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            // dsmc/generalBoundaryModels
            fileName dsmcGeneralBoundaryModelsTimePath
            (
                dsmcTimePath/"generalBoundaryModels"
            );

            if (!isDir(dsmcGeneralBoundaryModelsTimePath))
            {
                mkDir(dsmcGeneralBoundaryModelsTimePath);
            }

            forAll(generalBoundaryModels_, g)
            {
                if
                (
                    generalBoundaryModels_[g]->writeInTimeDir()
                 || generalBoundaryModels_[g]->writeInCase()
                )
                {
                    // directory: case/<timeDir>/uniform/controllers
                    // dsmc/generalBoundaryModels/<generalBoundaryModel>
                    fileName gMTimePath
                    (
                        dsmcGeneralBoundaryModelsTimePath
                        /gMFixedPathNames_[g]
                    );

                    if (!isDir(gMTimePath))
                    {
                        mkDir(gMTimePath);
                    }

                    // Creating directory for different zones but
                    // of the same model
                    const word& patchName =
                        generalBoundaryModels_[g]->patchName();

                    // directory: case/<timeDir>/uniform/controllers
                    // dsmc/generalBoundaryModels
                    // <generalBoundaryModel>/<patchName>
                    fileName patchTimePath(gMTimePath/patchName);

                    if (!isDir(patchTimePath))
                    {
                        mkDir(patchTimePath);
                    }

                    timePathNames[g] = patchTimePath;

                    generalBoundaryModels_[g]->output
                    (
                        gMFixedPathNames_[g], timePathNames[g]
                    );
                }
            }
        }
    }

    {
        patchBoundaryList_.clear();

        patchBoundaryList_ =
            dsmcBoundariesDict_.lookup("dsmcPatchBoundaries");

        forAll(patchBoundaryModels_, p)
        {
            const entry& boundaryI = patchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();

            patchBoundaryModels_[p]->updateProperties(boundaryIDict);
        }
    }

    {
        cyclicBoundaryList_.clear();

        cyclicBoundaryList_ =
            dsmcBoundariesDict_.lookup("dsmcCyclicBoundaries");

        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();

            cyclicBoundaryModels_[c]->updateProperties(boundaryIDict);
        }
    }

    {
        generalBoundaryList_.clear();

        generalBoundaryList_ =
            dsmcBoundariesDict_.lookup("dsmcGeneralBoundaries");

        forAll(generalBoundaryModels_, g)
        {
            const entry& boundaryI = generalBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();

            generalBoundaryModels_[g]->updateProperties(boundaryIDict);
        }
    }
}


label dsmcBoundaries::nPatchBoundaryModels() const
{
    return nPatchBoundaryModels_;
}


label dsmcBoundaries::nCyclicBoundaryModels() const
{
    return nCyclicBoundaryModels_;
}


label dsmcBoundaries::nGeneralBoundaryModels() const
{
    return nGeneralBoundaryModels_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
