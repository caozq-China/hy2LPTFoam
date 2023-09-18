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

#include "solidBoundaries.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::solidBoundaries::dictName("solidBoundariesDict");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Null Constructor 
solidBoundaries::solidBoundaries
(
    const Time& t,
    const polyMesh& mesh
)
:    
    time_(t),
    solidBoundariesDict_
    (
        IOobject
        (
            solidBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    nSolidPatchBoundaryModels_(0),
    nSolidCyclicBoundaryModels_(0),
    nSolidGeneralBoundaryModels_(0),

    solidPatchBoundaryList_(),
    solidPatchBoundaryNames_(),
    solidPatchBoundaryIds_(),
    solidPBFixedPathNames_(),
    solidPatchBoundaryModels_(),
    solidPatchToModelId_(mesh.boundaryMesh().size(), -1),

    solidCyclicBoundaryList_(),
    solidCyclicBoundaryNames_(),
    solidCyclicBoundaryIds_(),
    solidCMFixedPathNames_(),
    solidCyclicBoundaryModels_(),
    solidCyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    solidGeneralBoundaryList_(),
    solidGeneralBoundaryNames_(),
    solidGeneralBoundaryIds_(),
    solidGMFixedPathNames_(),
    solidGeneralBoundaryModels_()
{}


//- Constructor for dsmcFoam
solidBoundaries::solidBoundaries
(
    const Time& t,
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc
)
:
    time_(t),
    solidBoundariesDict_
    (
        IOobject
        (
            solidBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nSolidPatchBoundaryModels_(0),
    nSolidCyclicBoundaryModels_(0),
    nSolidGeneralBoundaryModels_(0),


    solidPatchBoundaryList_(solidBoundariesDict_.lookup("solidPatchBoundaries")),
    solidPatchBoundaryNames_(solidPatchBoundaryList_.size()),
    solidPatchBoundaryIds_(solidPatchBoundaryList_.size()),
    solidPBFixedPathNames_(solidPatchBoundaryList_.size()),
    solidPatchBoundaryModels_(solidPatchBoundaryList_.size()),
    solidPatchToModelId_(mesh.boundaryMesh().size(), -1),

    solidCyclicBoundaryList_(solidBoundariesDict_.lookup("solidCyclicBoundaries")),
    solidCyclicBoundaryNames_(solidCyclicBoundaryList_.size()),
    solidCyclicBoundaryIds_(solidCyclicBoundaryList_.size()),
    solidCMFixedPathNames_(solidCyclicBoundaryList_.size()),
    solidCyclicBoundaryModels_(solidCyclicBoundaryList_.size()),
    solidCyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),
    
    solidGeneralBoundaryList_(solidBoundariesDict_.lookup("solidGeneralBoundaries")),
    solidGeneralBoundaryNames_(solidGeneralBoundaryList_.size()),
    solidGeneralBoundaryIds_(solidGeneralBoundaryList_.size()),
    solidGMFixedPathNames_(solidGeneralBoundaryList_.size()),
    solidGeneralBoundaryModels_(solidGeneralBoundaryList_.size())
{
    Info << "Creating the solid boundary models: " << nl << endl;

    //- patch boundaries

    if( solidPatchBoundaryModels_.size() > 0 )
    {
        forAll(solidPatchBoundaryModels_, p)
        {
            const entry& boundaryI = solidPatchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            solidPatchBoundaryModels_[p] = autoPtr<solidPatchBoundary>
            (
                solidPatchBoundary::New(mesh, spc, boundaryIDict)
            );
    
            solidPatchBoundaryNames_[p] = solidPatchBoundaryModels_[p]->type();
            solidPatchBoundaryIds_[p] = p;
            ++nSolidPatchBoundaryModels_;
        }
    }

    checkSolidPatchBoundaryModels(mesh);

    //- cyclic boundaries

    if(solidCyclicBoundaryModels_.size() > 0 )
    {
        forAll(solidCyclicBoundaryModels_, c)
        {
            const entry& boundaryI = solidCyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            solidCyclicBoundaryModels_[c] = autoPtr<solidCyclicBoundary>
            (
                solidCyclicBoundary::New(mesh, spc, boundaryIDict)
            );
    
            solidCyclicBoundaryNames_[c] = solidCyclicBoundaryModels_[c]->type();
            solidCyclicBoundaryIds_[c] = c;
            ++nSolidCyclicBoundaryModels_;
        }
    }

    checkSolidCyclicBoundaryModels(mesh);
    
    //- solid general boundaries
 
    if(solidGeneralBoundaryModels_.size() > 0 )
    {
        forAll(solidGeneralBoundaryModels_, g)
        {
            const entry& boundaryI = solidGeneralBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            solidGeneralBoundaryModels_[g] = autoPtr<solidGeneralBoundary>
            (
                solidGeneralBoundary::New(mesh, spc, boundaryIDict)
            );
    
            solidGeneralBoundaryNames_[g] = solidGeneralBoundaryModels_[g]->type();
            solidGeneralBoundaryIds_[g] = g;
            ++nSolidGeneralBoundaryModels_;
        }
    }

    //- creating directories
    if(nSolidPatchBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/solidPhase
        fileName solidBoundariesPath(boundariesPath/"solid");

        if( !isDir(solidBoundariesPath) )
        {
            mkDir(solidBoundariesPath);
        }

        // directory: case/boundaries/solidPhase/patchBoundaryModels
        fileName patchBoundaryModelsPath(solidBoundariesPath/"patchBoundaryModels");
    
        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);    
        }

        forAll(solidPatchBoundaryModels_, p)
        {
            if(solidPatchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/dsmc/patchBoundaryModels/<patchBoundaryModel>
                fileName solidPatchBoundaryModelPath(patchBoundaryModelsPath/solidPatchBoundaryNames_[p]);

                if (!isDir(solidPatchBoundaryModelPath))
                {
                    mkDir(solidPatchBoundaryModelPath);    
                }
    
                const word& patchName = solidPatchBoundaryModels_[p]->patchName();

                // directory: case/controllers/dsmc/patchBoundaryModels/<patchBoundaryModel>/<patchName>    
                fileName patchPath(solidPatchBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                solidPBFixedPathNames_[p] = patchPath;
            }
        }
    }


    //- creating directories
    if(nSolidCyclicBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/dsmc
        fileName solidBoundariesPath(boundariesPath/"solid");

        if( !isDir(solidBoundariesPath) )
        {
            mkDir(solidBoundariesPath);
        }

        // directory: case/boundaries/dsmc/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath(solidBoundariesPath/"cyclicBoundaryModels");
    
        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);    
        }

        forAll(solidCyclicBoundaryModels_, c)
        {
            if(solidCyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/dsmc/cyclicBoundaryModels/<cyclicBoundaryModel>
                fileName solidCyclicBoundaryModelPath(cyclicBoundaryModelsPath/solidCyclicBoundaryNames_[c]);

                if (!isDir(solidCyclicBoundaryModelPath))
                {
                    mkDir(solidCyclicBoundaryModelPath);    
                }
    
                const word& patchName = solidCyclicBoundaryModels_[c]->patchName();

                // directory: case/controllers/dsmc/cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>      
                fileName patchPath(solidCyclicBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                solidCMFixedPathNames_[c] = patchPath;
            }
        }
    }
    
    //- creating solid directories
    if(nSolidGeneralBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/solidPhase
        fileName solidBoundariesPath(boundariesPath/"solid");

        if( !isDir(solidBoundariesPath) )
        {
            mkDir(solidBoundariesPath);
        }

        // directory: case/boundaries/solidPhase/cyclicBoundaryModels
        fileName solidGeneralBoundaryModelsPath(solidBoundariesPath/"solidGeneralBoundaryModels");
    
        if (!isDir(solidGeneralBoundaryModelsPath))
        {
            mkDir(solidGeneralBoundaryModelsPath);    
        }

        forAll(solidGeneralBoundaryModels_, g)
        {
            if(solidGeneralBoundaryModels_[g]->writeInCase())
            {
                // directory: case/boundaries/dsmc/generalBoundaryModels/<generalBoundaryModel>
                fileName solidGeneralBoundaryModelPath(solidGeneralBoundaryModelsPath/solidGeneralBoundaryNames_[g]);

                if (!isDir(solidGeneralBoundaryModelPath))
                {
                    mkDir(solidGeneralBoundaryModelPath);    
                }
    
                const word& patchName = solidGeneralBoundaryModels_[g]->patchName();

                // directory: case/controllers/dsmc/generalBoundaryModels/<generalBoundaryModel>/<patchName>      
                fileName patchPath(solidGeneralBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                solidGMFixedPathNames_[g] = patchPath;
            }
        }
    }

}


void solidBoundaries::checkSolidCyclicBoundaryModels(const polyMesh& mesh)
{
    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
    
        if(isA<cyclicPolyPatch>(patch))
        {
            label patchIndex = patch.index();

            forAll(solidCyclicBoundaryModels_, c)
            {
                const label& patchId = solidCyclicBoundaryModels_[c]->patchId();
 
                if(patchIndex == patchId)
                {
                    ++nPolyPatches;
                    solidCyclicBoundaryToModelId_[patchi] = c;
                }
            }
        }
    }

    if(nPolyPatches != nSolidCyclicBoundaryModels_)
    {
        FatalErrorInFunction
            << nl
            << " Number of cyclic boundary models = "  << nSolidCyclicBoundaryModels_ 
            << " chosen in the solidBoundariesDict are inconsistent." 
            << abort(FatalError);
    }
}


void solidBoundaries::checkSolidPatchBoundaryModels(const polyMesh& mesh)
{
    //- check that all poly-patches defined within blockMeshDict,
    //  each have one model.

    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
    
        if (!polyPatch::constraintType(patch.type()))
        {
            ++nPolyPatches;

            label patchIndex = patch.index();

            label nPatches = 0;

            forAll(solidPatchBoundaryModels_, p)
            {
                const label patchId = solidPatchBoundaryModels_[p]->patchId();
 
                if(patchIndex == patchId)
                {
                    ++nPatches;
                    solidPatchToModelId_[patchi] = p;
                }
            }

            if(nPatches > 1)
            {
                FatalIOErrorInFunction(solidBoundariesDict_)
                    << nl
                    << " Only one patch boundary model per poly-patch, [name: "
                    << patch.name()
                    << "]. No of models chosen for this patch are: " 
                    << nPatches  << ", in " 
                    << mesh.time().system()/"solidBoundariesDict"
                    << abort(FatalError);
            }
        }
    }

//     Pout << "patchToModelId_: " << patchToModelId_ << endl;

    if(nPolyPatches != nSolidPatchBoundaryModels_)
    {
        FatalIOErrorInFunction(solidBoundariesDict_)
            << nl
            << " Number of poly-patches = "  << nPolyPatches 
            << " are not equal to the number of patch "
            << "models = " << nSolidPatchBoundaryModels_
            << abort(FatalError);
    }
}


void solidBoundaries::updateTimeInfo()
{
    for (auto& model : solidGeneralBoundaryModels_)
    {
        model->updateTime();
    }
}

void solidBoundaries::setInitialConfig()
{
    for (auto& model : solidPatchBoundaryModels_)
    {
        model->initialConfiguration();
    }

    for (auto& model : solidCyclicBoundaryModels_)
    {
        model->initialConfiguration();
    }

    for (auto& model : solidGeneralBoundaryModels_)
    {
        model->initialConfiguration();
    }
}

void solidBoundaries::calculateProps()
{
    for (auto& model : solidPatchBoundaryModels_)
    {
        model->calculateProperties();
    }

    for (auto& model : solidCyclicBoundaryModels_)
    {
        model->calculateProperties();
    }
    
    for (auto& model : solidGeneralBoundaryModels_)
    {
        model->calculateProperties();
    }
}

void solidBoundaries::releaseParticlesFromWall()
{
    for (auto& model : solidPatchBoundaryModels_)
    {
        model->injectParticlesFromWall();
    }

}

// impose model after calculation of forces
void solidBoundaries::controlBeforeMove()
{
    for (auto& model : solidGeneralBoundaryModels_)
    {
        model->controlParcelsBeforeMove();
    }
}


void solidBoundaries::controlBeforeCollisions()
{
    for (auto& model : solidGeneralBoundaryModels_)
    {
        model->controlParcelsBeforeCollisions();
    }
}

void solidBoundaries::controlAfterCollisions()
{
    for (auto& model : solidGeneralBoundaryModels_)
    {
        model->controlParcelsAfterCollisions();
    }
}


//- output
void solidBoundaries::outputResults()
{
    if (!time_.writeTime())
    {
        return;
    }

//     if(runTime.outputTime())
//     {
        {
            List<fileName> timePathNames(solidPBFixedPathNames_.size()); 
    
            if(nSolidPatchBoundaryModels_ > 0)
            {
                // directory: case/<timeDir>/uniform
                fileName uniformTimePath(time_.path()/time_.timeName()/"uniform");
            
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

                // directory: case/<timeDir>/uniform/boundaries/solid
                fileName solidTimePath(boundariesTimePath/"solid");
            
                if (!isDir(solidTimePath))
                {
                    mkDir(solidTimePath);    
                }

                // directory: case/<timeDir>/uniform/boundaries/solidPhase/patchBoundaryModels
                fileName solidPatchBoundaryModelsTimePath(solidTimePath/"solidPatchBoundaryModels");
            
                if (!isDir(solidPatchBoundaryModelsTimePath))
                {
                    mkDir(solidPatchBoundaryModelsTimePath);    
                }

                forAll(solidPatchBoundaryModels_, p)
                {
                    if
                    (
                        solidPatchBoundaryModels_[p]->writeInTimeDir() ||
                        solidPatchBoundaryModels_[p]->writeInCase()
                    )
                    {
						
//                         // directory: case/<timeDir>/uniform/controllers/dsmc/patchBoundaryModels/<patchBoundaryModel>
//                         fileName pBTimePath(dsmcPatchBoundaryModelsTimePath/pBFixedPathNames_[p]);
// 
//                         if(!isDir(pBTimePath))
//                         {
//                             mkDir(pBTimePath);
//                         }

                        //- creating directory for different zones but of the same model
                        const word& patchName = solidPatchBoundaryModels_[p]->patchName();

                        // directory: case/<timeDir>/uniform/controllers/dsmc/patchBoundaryModels/<patchBoundaryModel>/<patchName>
                        fileName patchTimePath(solidPatchBoundaryModelsTimePath/patchName);

                        if (!isDir(patchTimePath))
                        {
                            mkDir(patchTimePath);
                        }

                        timePathNames[p] = patchTimePath;

                        solidPatchBoundaryModels_[p]->output(solidPBFixedPathNames_[p], timePathNames[p]);
// 						patchBoundaryModels_[p]->writeField();
                    }
                }
            }
        } 

        {
            List<fileName> timePathNames(solidCMFixedPathNames_.size()); 
    
            if(nSolidCyclicBoundaryModels_ > 0)
            {
                // directory: case/<timeDir>/uniform
                fileName uniformTimePath(time_.path()/time_.timeName()/"uniform");
            
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

                // directory: case/<timeDir>/uniform/boundaries/solid
                fileName solidTimePath(boundariesTimePath/"solid");
            
                if (!isDir(solidTimePath))
                {
                    mkDir(solidTimePath);    
                }

                // directory: case/<timeDir>/uniform/boundaries/dsmc/cyclicBoundaryModels
                fileName solidCyclicBoundaryModelsTimePath(solidTimePath/"solidCyclicBoundaryModels");
            
                if (!isDir(solidCyclicBoundaryModelsTimePath))
                {
                    mkDir(solidCyclicBoundaryModelsTimePath);    
                }

                forAll(solidCyclicBoundaryModels_, c)
                {
                    if
                    (
                        solidCyclicBoundaryModels_[c]->writeInTimeDir() ||
                        solidCyclicBoundaryModels_[c]->writeInCase()
                    )
                    {
                        // directory: case/<timeDir>/uniform/controllers/dsmc/cyclicBoundaryModels/<cyclicBoundaryModel>
                        fileName cMTimePath(solidCyclicBoundaryModelsTimePath/solidCMFixedPathNames_[c]);

                        if(!isDir(cMTimePath))
                        {
                            mkDir(cMTimePath);
                        }

                        //- creating directory for different zones but of the same model
                        const word& patchName = solidCyclicBoundaryModels_[c]->patchName();

                        // directory: case/<timeDir>/uniform/controllers/dsmc/cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>
                        fileName patchTimePath(cMTimePath/patchName);

                        if (!isDir(patchTimePath))
                        {
                            mkDir(patchTimePath);
                        }

                        timePathNames[c] = patchTimePath;

                        solidCyclicBoundaryModels_[c]->output(solidCMFixedPathNames_[c], timePathNames[c]);
                    }
                }
            }
        } 

        
        {
            List<fileName> timePathNames(solidGMFixedPathNames_.size()); 
    
            if(nSolidGeneralBoundaryModels_ > 0)
            {
                // directory: case/<timeDir>/uniform
                fileName uniformTimePath(time_.path()/time_.timeName()/"uniform");
            
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

                // directory: case/<timeDir>/uniform/boundaries/solid
                fileName solidTimePath(boundariesTimePath/"solid");
            
                if (!isDir(solidTimePath))
                {
                    mkDir(solidTimePath);    
                }

                // directory: case/<timeDir>/uniform/boundaries/solidPhase/generalBoundaryModels
                fileName solidGeneralBoundaryModelsTimePath(solidTimePath/"solidGeneralBoundaryModels");
            
                if (!isDir(solidGeneralBoundaryModelsTimePath))
                {
                    mkDir(solidGeneralBoundaryModelsTimePath);    
                }

                forAll(solidGeneralBoundaryModels_, g)
                {
                    if
                    (
                        solidGeneralBoundaryModels_[g]->writeInTimeDir() ||
                        solidGeneralBoundaryModels_[g]->writeInCase()
                    )
                    {
                        // directory: case/<timeDir>/uniform/controllers/dsmc/generalBoundaryModels/<generalBoundaryModel>
                        fileName solidGMTimePath(solidGeneralBoundaryModelsTimePath/solidGMFixedPathNames_[g]);

                        if(!isDir(solidGMTimePath))
                        {
                            mkDir(solidGMTimePath);
                        }

                        //- creating directory for different zones but of the same model
                        const word& patchName = solidGeneralBoundaryModels_[g]->patchName();

                        // directory: case/<timeDir>/uniform/controllers/dsmc/generalBoundaryModels/<generalBoundaryModel>/<patchName>
                        fileName patchTimePath(solidGMTimePath/patchName);

                        if (!isDir(patchTimePath))
                        {
                            mkDir(patchTimePath);
                        }

                        timePathNames[g] = patchTimePath;

                        solidGeneralBoundaryModels_[g]->output(solidGMFixedPathNames_[g], timePathNames[g]);
                    }
                }
            }
        }

        // RE-READ DICTIONARIES FOR MODIFIED PROPERTIES (RUN-TIME SELECTION)

        {
            solidPatchBoundaryList_.clear();
        
            solidPatchBoundaryList_ = solidBoundariesDict_.lookup("solidPatchBoundaries");
        
            forAll(solidPatchBoundaryModels_, p)
            {
                const entry& boundaryI = solidPatchBoundaryList_[p];
                const dictionary& boundaryIDict = boundaryI.dict();
    
                solidPatchBoundaryModels_[p]->updateProperties(boundaryIDict);
            }
        }

        {
            solidCyclicBoundaryList_.clear();
        
            solidCyclicBoundaryList_ = solidBoundariesDict_.lookup("solidCyclicBoundaries");
        
            forAll(solidCyclicBoundaryModels_, c)
            {
                const entry& boundaryI = solidCyclicBoundaryList_[c];
                const dictionary& boundaryIDict = boundaryI.dict();
    
                solidCyclicBoundaryModels_[c]->updateProperties(boundaryIDict);
            }
        }


        
        {
            solidGeneralBoundaryList_.clear();
        
            solidGeneralBoundaryList_ = solidBoundariesDict_.lookup("solidGeneralBoundaries");
        
            forAll(solidGeneralBoundaryModels_, g)
            {
                const entry& boundaryI = solidGeneralBoundaryList_[g];
                const dictionary& boundaryIDict = boundaryI.dict();
    
                solidGeneralBoundaryModels_[g]->updateProperties(boundaryIDict);
            }
        }
//     }
}

label solidBoundaries::nSolidPatchBoundaryModels() const
{
    return nSolidPatchBoundaryModels_;
}

label solidBoundaries::nSolidCyclicBoundaryModels() const
{
    return nSolidCyclicBoundaryModels_;
}

label solidBoundaries::nSolidGeneralBoundaryModels() const
{
    return nSolidGeneralBoundaryModels_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
