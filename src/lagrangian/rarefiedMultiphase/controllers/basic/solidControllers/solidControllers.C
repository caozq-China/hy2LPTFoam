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

#include "solidControllers.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor 
solidControllers::solidControllers
(
    const Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    solidControllersDict_
    (
        IOobject
        (
            "solidControllersDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
//     nFluxControllers_(0),
    stateControllersList_(),
    sCNames_(),
    sCIds_(),
    sCFixedPathNames_(),
    stateControllers_()

//     fluxControllersList_(),
//     fCNames_(),
//     fCIds_(),
//     fCFixedPathNames_(),
//     fluxControllers_()
{}

//- Constructor for mdFOAM
solidControllers::solidControllers
(
    const Time& t,
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud
)
:
    time_(t),
    solidControllersDict_
    (
        IOobject
        (
            "solidControllersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
//     nFluxControllers_(0),
    stateControllersList_(solidControllersDict_.lookup("solidStateControllers")),
    sCNames_(stateControllersList_.size()),
    sCIds_(stateControllersList_.size()),
    sCFixedPathNames_(stateControllersList_.size()),
    stateControllers_(stateControllersList_.size())

//     fluxControllersList_(solidControllersDict_.lookup("solidFluxControllers")),
//     fCNames_(fluxControllersList_.size()),
//     fCIds_(fluxControllersList_.size()),
//     fCFixedPathNames_(fluxControllersList_.size()),
//     fluxControllers_(fluxControllersList_.size())
{

    Info << nl << "Creating solidControllers" << nl << endl;

    //- state solidControllers
    nStateControllers_ = stateControllers_.size();

    forAll(stateControllers_, sC)
    {
        const entry& solidControllersI = stateControllersList_[sC];
        const dictionary& solidControllersIDict = solidControllersI.dict();

        stateControllers_[sC] = solidStateController::New(cloud.mesh(), cloud, solidControllersIDict);

        sCNames_[sC] = stateControllers_[sC]->type();
        sCIds_[sC] = sC;
    }

//     //- flux solidControllers
// 
//     if(fluxControllers_.size() > 0 )
//     {
//         forAll(fluxControllers_, fC)
//         {
//             const entry& solidControllersI = fluxControllersList_[fC];
//     
//             const dictionary& solidControllersIDict = solidControllersI.dict();
//     
//             fluxControllers_[fC] = autoPtr<solidFluxController>
//             (
//                 solidFluxController::New(time_, cloud, solidControllersIDict)
//             );
//     
//             fCNames_[fC] = fluxControllers_[fC]->type();
//             fCIds_[fC] = fC;
//     
//             nFluxControllers_++;
//         }
//     }

    // creating directories for state controllers
    if(nStateControllers_ > 0)
    {
        // directory: case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if( !isDir(controllersPath) )
        {
            mkDir(controllersPath);
        }

        // directory: case/controllers/solid
        fileName solidControllersPath(controllersPath/"solid");

        if( !isDir(solidControllersPath) )
        {
            mkDir(solidControllersPath);
        }

        // directory: case/controllers/solid/stateControllers
        fileName stateControllersPath(solidControllersPath/"stateControllers");
    
        if (!isDir(stateControllersPath))
        {
            mkDir(stateControllersPath);    
        }

        forAll(stateControllers_, sC)
        {
            if(stateControllers_[sC]->writeInCase())
            {
                // directory: case/controllers/solid/stateControllers/<stateControllerModel>
                fileName stateControllerPath(stateControllersPath/sCNames_[sC]);

                if (!isDir(stateControllerPath))
                {
                    mkDir(stateControllerPath);    
                }
    
                const word& regionName = stateControllers_[sC]->regionName();

                // directory: case/controllers/solid/stateControllers/<stateControllerModel>/<cellZoneName>    
                fileName zonePath(stateControllerPath/regionName);
   
                if (!isDir(zonePath))
                {
                    mkDir(zonePath);    
                }
    
                sCFixedPathNames_[sC] = zonePath;
            }
        }
    }

    // creating directories for flux controllers
//     if(nFluxControllers_ > 0)
//     {
//         // directory: case/controllers
//         fileName controllersPath(time_.path()/"controllers");
// 
//         if( !isDir(controllersPath) )
//         {
//             mkDir(controllersPath);
//         }
// 
//         // directory: case/controllers/solid
//         fileName solidControllersPath(time_.path()/"solid");
// 
//         if( !isDir(solidControllersPath) )
//         {
//             mkDir(solidControllersPath);
//         }
// 
//         // directory: case/controllers/solid/fluxControllers
//         fileName fluxControllersPath(solidControllersPath/"fluxControllers");
//     
//         if (!isDir(fluxControllersPath))
//         {
//             mkDir(fluxControllersPath);    
//         }
// 
//         forAll(fluxControllers_, fC)
//         {
//             if(fluxControllers_[fC]->writeInCase())
//             {
//                 // directory: case/controllers/solid/fluxControllers/<fluxControllerModel>
//                 fileName fluxControllerPath(fluxControllersPath/fCNames_[fC]);
//     
//                 if (!isDir(fluxControllerPath))
//                 {
//                     mkDir(fluxControllerPath);    
//                 }
// 
//                 const word& regionName = fluxControllers_[fC]->regionName();
//     
//                 // directory: case/controllers/solid/fluxControllers/<fluxControllerModel>/<faceZoneName>
//                 fileName zonePath(fluxControllerPath/regionName);
//     
//                 if (!isDir(zonePath))
//                 {
//                     mkDir(zonePath);    
//                 }
// 
//                 fCFixedPathNames_[fC] = zonePath;
//             }
//         }
//     }
}


//- initial configuration
//- call this function after the solidMoleculeCloud is completely initialised
void solidControllers::initialConfig()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->initialConfiguration();
    }

//     forAll(fluxControllers_, fC)
//     {
//         fluxControllers_[fC]->initialConfiguration();
//     }
}

        //- different control stages 
void solidControllers::controlBeforeMove()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlSolidParticlesBeforeMove();
    }
}

void solidControllers::controlAfterMove()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlSolidParticlesAfterMove();
    }
}


void solidControllers::controlBeforeCollisions()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlSolidParticlesBeforeCollisions();
    }
}

void solidControllers::controlAfterCollisions()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlSolidParticlesAfterCollisions();
    }
}



//- calculate properties -- call this at the end of the MD time-step.
void solidControllers::calculateProps()
{
    forAll(stateControllers_, sC)
    {
//         Info << "error: " << sCNames_[sC] << endl;
        stateControllers_[sC]->calculateProperties();
    }

//     forAll(fluxControllers_, fC)
//     {
// //         Info << "error: " << sCNames_[sC] << endl;
//         fluxControllers_[fC]->calculateProperties();
//     }
}

//- this function is to be called at the beginning of the MD time-step. 
//  since we have placed a non-referenced time-data class in the state-controller class.
void solidControllers::updateTimeInfo()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->updateTime();
    }

//     forAll(fluxControllers_, fC)
//     {
//         fluxControllers_[fC]->updateTime();
//     }
}

//- output -- call this function at the end of the MD time-step
void solidControllers::outputResults() 
{
//     const Time& runTime = time_;
    if (!time_.writeTime())
    {
        return;
    }

//     if(runTime.outputTime())
//     {
        // -- creating a set of directories in the current time directory
        {
            List<fileName> timePathNames(sCFixedPathNames_.size()); 
    
            if(nStateControllers_ > 0)
            {
                // directory: case/<timeDir>/uniform
                fileName uniformTimePath(time_.path()/time_.timeName()/"uniform");
            
                if (!isDir(uniformTimePath))
                {
                    mkDir(uniformTimePath);
                }
                
                // directory: case/<timeDir>/uniform/controllers
                fileName controllersTimePath(uniformTimePath/"controllers");

                if (!isDir(controllersTimePath))
                {
                    mkDir(controllersTimePath);
                }

                // directory: case/<timeDir>/uniform/controllers/dsmc
                fileName solidTimePath(controllersTimePath/"solid");

                if (!isDir(solidTimePath))
                {
                    mkDir(solidTimePath);
                }

                // directory: case/<timeDir>/uniform/controllers/dsmc/
                fileName solidStateControllersTimePath
                (
                    solidTimePath/"stateControllers"
                );

                if (!isDir(solidStateControllersTimePath))
                {
                    mkDir(solidStateControllersTimePath);
                }
                
                forAll(stateControllers_, sC)
                {
                    if (stateControllers_[sC]->writeInTimeDir())
                    {
                        // directory:
                        // case/<timeDir>/uniform/controllers/
                        // dsmc/<stateControllerModel>
                        fileName sCTimePath =
                            solidStateControllersTimePath/sCNames_[sC];

                        if (!isDir(sCTimePath))
                        {
                            mkDir(sCTimePath);
                        }

                        // creating directory for different zones but
                        // of the same model
                        const word& regionName =
                            stateControllers_[sC]->regionName();

                        // directory: case/<timeDir>/uniform/controllers/
                        // dsmc/<stateControllerModel>/
                        // <cellZoneName>
                        fileName zoneTimePath(sCTimePath/regionName);

                        if (!isDir(zoneTimePath))
                        {
                            mkDir(zoneTimePath);
                        }

                        timePathNames[sC] = zoneTimePath;
                    }
                }
            }
        
            // -- write out data (do not comment this out)
            forAll(stateControllers_, sC)
            {
                stateControllers_[sC]->output(sCFixedPathNames_[sC], timePathNames[sC]);
            }
        }

//         {
//             List<fileName> timePathNames(fCFixedPathNames_.size());
//     
//             if(nFluxControllers_ > 0)
//             {
//                 // directory: case/<timeDir>/uniform
//                 fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
//             
//                 if (!isDir(uniformTimePath))
//                 {
//                     mkDir(uniformTimePath);
//                 }
//     
//                 if(fluxControllers_.size() > 0)
//                 {
//                    // directory: case/<timeDir>/uniform/controllers
//                     fileName controllersTimePath(uniformTimePath/"controllers");
// 
//                     if (!isDir(controllersTimePath))
//                     {
//                         mkDir(controllersTimePath);
//                     }
// 
//                     // directory: case/<timeDir>/uniform/controllers/solid
//                     fileName solidTimePath(controllersTimePath/"solid");
//                 
//                     if (!isDir(solidTimePath))
//                     {
//                         mkDir(solidTimePath);    
//                     }
// 
//                     // directory: case/<timeDir>/uniform/fluxControllers
//                     fileName solidControllersTimePath(solidTimePath/"fluxControllers");
//                 
//                     if (!isDir(solidControllersTimePath))
//                     {
//                         mkDir(solidControllersTimePath);    
//                     }
// 
//                     forAll(fluxControllers_, fC)
//                     {
//                         if
//                         (
//                             fluxControllers_[fC]->writeInTimeDir()
//                         )
//                         {
//                             // directory: case/<timeDir>/uniform/controllers/solid/<fluxControllerModel>
//                             fileName fCTimePath(solidControllersTimePath/fCNames_[fC]);
//         
//                             if(!isDir(fCTimePath))
//                             {
//                                 mkDir(fCTimePath);
//                             }
//         
//                             const word& regionName = fluxControllers_[fC]->regionName();
// 
//                             // directory: case/<timeDir>/uniform/controllers/solid/<fluxControllerModel>  <faceZoneName>      
//                             fileName zoneTimePath(fCTimePath/regionName);
//             
//                             if (!isDir(zoneTimePath))
//                             {
//                                 mkDir(zoneTimePath);    
//                             }
// 
//                             timePathNames[fC] = zoneTimePath;
//                         }
//                     }
//                 }
//             }
// 
//             // -- write out data (do not comment this out)
//             forAll(fluxControllers_, fC)
//             {
//                 fluxControllers_[fC]->output(fCFixedPathNames_[fC], timePathNames[fC]);
//             }
//         }

        // RE-READ DICTIONARIES FOR MODIFIED PROPERTIES (RUN-TIME SELECTION)

        {
            stateControllersList_.clear();
        
            stateControllersList_ = solidControllersDict_.lookup("solidStateControllers");
        
            forAll(stateControllers_, sC)
            {
                const entry& solidControllersI = stateControllersList_[sC];
                const dictionary& solidControllersIDict = solidControllersI.dict();
    
                stateControllers_[sC]->updateProperties(solidControllersIDict);
            }
        }

//         {
//             fluxControllersList_.clear();
//         
//             fluxControllersList_ = solidControllersDict_.lookup("solidFluxControllers");
//         
//             forAll(fluxControllers_, fC)
//             {
//                 const entry& solidControllersI = fluxControllersList_[fC];
//                 const dictionary& solidControllersIDict = solidControllersI.dict();
//     
//                 fluxControllers_[fC]->updateProperties(solidControllersIDict);
//             }
//         }
//     }
}


label solidControllers::nStateControllers() const
{
    return nStateControllers_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
