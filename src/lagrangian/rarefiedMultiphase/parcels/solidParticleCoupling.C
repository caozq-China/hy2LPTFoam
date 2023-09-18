/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "solidParticleCoupling.H"
#include "solidParticleCouplingCloud.H"
#include "meshTools.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTemplateTypeNameAndDebug(Cloud<solidParticleCoupling>, 0);
// }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::solidParticleCoupling::move
(
    solidParticleCouplingCloud& spc,
    trackingData& td,
    const scalar trackTime
)
{
//     trackingData::solidParticleCouplingCloud::solidParticleCoupling& p = static_cast<trackingData::solidParticleCouplingCloud::solidParticleCoupling&>(*this);
    
    switch (td.part())
    {
        case trackingData::LinearTrack:
        {
            td.switchProcessor = false;
            td.keepParticle = true;
            
            const polyMesh& mesh = spc.pMesh();
//             const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

            if(newSolidParticle() == 1)
            {
                label& nP = newSolidParticle();
                Random& rndGenS = spc.rndGenS();
                stepFraction() = rndGenS.sample01<scalar>(); 
                nP = 0;
            }
            
//             scalar tEnd = (1.0 - stepFraction())*trackTime;
//             const scalar dtMax = tEnd;

            vector Utracking = U_;
        
            while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
            {
                Utracking = U_;
                
                // Apply correction to velocity to constrain tracking for
                // reduced - D cases
                meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
                
                // Deviation from the mesh centre for reduced - D cases
                const vector d = deviationFromMeshCentre();
                
                const scalar f = 1 - stepFraction();
                trackToAndHitFace(f*trackTime*Utracking - d, f, spc, td);
                
//                 if (onFace())
//                 {
//                     //-monitoring flux properties
//                     spc.tracker().updateFields(*this);
//                     spc.functions().postFace(*this, td.keepParticle);
//                 }

                if (onBoundaryFace() && td.keepParticle)
                {
                    forAll(spc.boundaries().solidCyclicBoundaryModels(), c)
                    {
                        const labelList& faces =
                            spc.boundaries().solidCyclicBoundaryModels()[c]->allFaces();

                        if (faces.find(this->face()) != -1)
                        {
                            spc.boundaries().solidCyclicBoundaryModels()[c]->controlMol
                            (
                                *this, td
                            );
                        }
                    }
                }
            }
            
//             Info<<"Leave move()"<<endl;
            
            break;
        }
//         case trackingData::DampingNoTrack:
//         {
//             UCorrect_ = spc.mppicDampingModels().velocityCorrection(*this,trackTime);
//             
//             td.keepParticle = true;
//             td.switchProcessor = false;
//             
//             break;
//         }
//         case trackingData::PackingNoTrack:
//         {
//             UCorrect_ = spc.mppicPackingModels().velocityUpdateAndCorrection(*this,trackTime);
//             
//             td.keepParticle = true;
//             td.switchProcessor = false;
//             
//             break;
//         }
        case trackingData::CorrectTrack:
        {   
            td.switchProcessor = false;
            td.keepParticle = true;
            if (UCorrect_.x()!=0.0 || UCorrect_.y()!=0.0 || UCorrect_.z()!=0.0)
            {
                vector U = U_;
        
                scalar f = stepFraction();
                
                U_ = (1.0 - f)*UCorrect_;
                    
            // - - - - - - - - - - Begin parcelType::move()  - - - - - - - - - -//
//                 td.switchProcessor = false;
//                 td.keepParticle = true;
                
                const polyMesh& mesh = spc.pMesh();
    //             const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

                if(newSolidParticle() == 1)
                {
                    label& nP = newSolidParticle();
                    Random& rndGenS = spc.rndGenS();
                    stepFraction() = rndGenS.sample01<scalar>(); 
                    nP = 0;
                }
                
    //             scalar tEnd = (1.0 - stepFraction())*trackTime;
    //             const scalar dtMax = tEnd;

                vector Utracking = U_;
                
                while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
                {
                    Utracking = U_;
                    
                    // Apply correction to velocity to constrain tracking for
                    // reduced - D cases
                    meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
                    
                    // Deviation from the mesh centre for reduced - D cases
                    const vector d = deviationFromMeshCentre();
                    
                    const scalar f = 1 - stepFraction();
                    trackToAndHitFace(f*trackTime*Utracking - d, f, spc, td);
                    
    //                 if (onFace())
    //                 {
    //                     //-monitoring flux properties
    //                     spc.tracker().updateFields(*this);
    //                     spc.functions().postFace(*this, td.keepParticle);
    //                 }

                    if (onBoundaryFace() && td.keepParticle)
                    {
                        forAll(spc.boundaries().solidCyclicBoundaryModels(), c)
                        {
                            const labelList& faces =
                                spc.boundaries().solidCyclicBoundaryModels()[c]->allFaces();

                            if (faces.find(this->face()) != -1)
                            {
                                spc.boundaries().solidCyclicBoundaryModels()[c]->controlMol
                                (
                                    *this, td
                                );
                            }
                        }
                    }
                }
            // - - - - - - - - - - End parcelType::move()  - - - - - - - - - -//

                U_ = U + (stepFraction() - f)*UCorrect_;
            }
            break;
        }
    }
    
    return td.keepParticle;
}

bool Foam::solidParticleCoupling::hitPatch
(
    solidParticleCouplingCloud& spc,
    trackingData& td
)
{
    return false;
}

void Foam::solidParticleCoupling::hitProcessorPatch
(
    solidParticleCouplingCloud& spc,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::solidParticleCoupling::hitWallPatch
(
    solidParticleCouplingCloud& spc,
    trackingData& td
)
{
    //-find which patch has been hit
//     label patchIndex = wpp.index();

    const label& patchModelId = spc.boundaries().
    solidPatchToModelIds()[patch()];

    // apply a boundary model when a molecule collides with this poly patch
    spc.boundaries().solidPatchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}


void Foam::solidParticleCoupling::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::solidParticleCoupling::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


// Foam::scalar Foam::solidParticleCoupling::wallImpactDistance(const vector&) const
// {
//     return 0.5*constProperties_.dSolid(); 
// }

// ************************************************************************* //
// #include "solidParticleCouplingIO.C"

