/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Time-driven Hard Sphere

Description

\*----------------------------------------------------------------------------*/

#include "coarseGrainedHardSphere.H"
#include "addToRunTimeSelectionTable.H"
#include "labelListIOList.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coarseGrainedHardSphere, 0);

addToRunTimeSelectionTable
(solidCollisionDetection, coarseGrainedHardSphere, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
coarseGrainedHardSphere::coarseGrainedHardSphere
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidCollisionDetection(mesh, spc, dict),
    dict_(dict.subDict("gridBasedNeighborSearchProperties")),
    infoCounter_(0),
    nCorrectionSteps_(dict_.get<label>("nCorrectionSteps")),
    k_(dict_.get<scalar>("K")),
    particlePtr_(nullptr)
//     pLabels_()
{
//     pLabels_.setSize(spc.size());
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool coarseGrainedHardSphere::active() const
{
    return true;
}


void coarseGrainedHardSphere::initialConfiguration()
{
    //- build neighbor collision candidates for each particle
    labelListIOList ngbCollisionCandidates
    (
        IOobject
        (
            "ngbParticles",
            mesh_.time().constant(),
            mesh_,
            IOobject::LAZY_READ,
            IOobject::AUTO_WRITE
        )
    );
    ngbCollisionCandidates.setSize(spc_.size());

    labelListIOList ngbCellsList
    (
        IOobject
        (
            "ngbCells",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar distance;
    scalar threshold;

    label count=0;
//     scalar collisionDiameter_p = 0;
//     scalar collisionDiameter_pj = 0;
    //- find the possible collision candinate
    forAllIter(solidParticleCouplingCloud, spc_, iter)
    {
        distance = 0.0;
        threshold = 0.0;

        solidParticleCoupling& p = iter();
//         solidParticleCoupling* p = &iter();

        forAll(ngbCellsList[p.cell()],ngbCellI)
        {
            //- ngbCellsList[p->cell()][ngbCellI] is the exact neighbour cell ID
            const DynamicList<solidParticleCoupling*>& ngbCellsolidParcels=spc_.cellOccupancy()[ngbCellsList[p.cell()][ngbCellI]];

//             collisionDiameter_p = 0;
            forAll(ngbCellsolidParcels,i)
            {
//                 collisionDiameter_pj = 0;
                const solidParticleCoupling& pj = *ngbCellsolidParcels[i];
                distance = mag(p.position()-pj.position());

                if(spc_.dsmcCloudReference()->axisymmetric())
                {
                    threshold = k_*(p.D()*pow(p.RWF(),1/3)+pj.D()*pow(pj.RWF(),1/3))/2;
//                     collisionDiameter_p = p.D()*pow((p.RWF()*spc_.nSolidParticles()),1/3);
//                     collisionDiameter_pj = pj.D()*pow((pj.RWF()*spc_.nSolidParticles()),1/3);
                }
                else
                {
                    threshold = k_*(p.D()+pj.D())/2;
//                     collisionDiameter_p = p.D()*pow(spc_.nSolidParticles(),1/3);
//                     collisionDiameter_pj = pj.D()*pow(spc_.nSolidParticles(),1/3);
                }

//                 threshold = k_*(collisionDiameter_p+collisionDiameter_pj)/2;



                if(distance<threshold)
                {
                    //- append particle numerical sequence:1,2,3,4...
                    ngbCollisionCandidates[count].append(pj.numSeq());
                }
            }
        }

        //- sort the particle number sequence
        //-i.e. [23,10,13,5,6,...]
        //-To [5,6,10,13,23,...]
        sort(ngbCollisionCandidates[count]);
        count++;
    }

    ngbCollisionCandidates.write();

    //- build particle Ptr list
    if(!particlePtr_)
    {
        particlePtr_.reset(new List<solidParticleCoupling*>(spc_.size()));
    }
    else if(particlePtr_().size() != spc_.size())
    {
        particlePtr_().setSize(spc_.size());
    }

    auto& particlePtr = particlePtr_();

    particlePtr_().clear();

    for (solidParticleCoupling& p: spc_)
    {
        //- save particle address
        particlePtr.append(&p);
    }
}


void coarseGrainedHardSphere::collide ()
{
    if (!spc_.solidBinaryCollision().active())
    {
        return;
    }

//     const List<DynamicList<solidParticleCoupling*> >& cellOccupancy =
//                                                 spc_.cellOccupancy();

    const polyMesh& mesh = spc_.dsmcCloudReference()->mesh();

    labelListIOList ngbParticles
    (
        IOobject
        (
            "ngbParticles",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label count = 0;
    scalar particleCenterDistance;
    scalar twoParticleDiameter;
//     scalar collisionDiameter_p = 0;
//     scalar collisionDiameter_pj = 0;
    //- binary collisions
    forAllIter(solidParticleCouplingCloud, spc_, iter)
    {
//         solidParticleCoupling* p = &iter();
        solidParticleCoupling& p = iter();

        particleCenterDistance = 0;
        twoParticleDiameter = 0;
//         collisionDiameter_p = 0;

        forAll(ngbParticles[count],j)
        {
            //loop over neighbour particles in neighbor ngbCells
            //- to perform binary collisions
            solidParticleCoupling& pj = *particlePtr_()[ngbParticles[count][j]];

//             collisionDiameter_pj = 0;

            particleCenterDistance = mag(p.position()-pj.position());

            if(spc_.dsmcCloudReference()->axisymmetric())
            {
                twoParticleDiameter = (p.D()*pow(p.RWF(),1/3)+pj.D()*pow(pj.RWF(),1/3))/2;
//                 collisionDiameter_p = p.D()*pow((p.RWF()*spc_.nSolidParticles()),1/3);
//                 collisionDiameter_pj = pj.D()*pow((pj.RWF()*spc_.nSolidParticles()),1/3);
            }
            else
            {
                twoParticleDiameter = (p.D()+pj.D())/2;
//                 collisionDiameter_p = p.D()*pow(spc_.nSolidParticles(),1/3);
//                 collisionDiameter_pj = pj.D()*pow(spc_.nSolidParticles(),1/3);
            }

//             twoParticleDiameter = (collisionDiameter_p+collisionDiameter_pj)/2;



            if(particleCenterDistance<twoParticleDiameter)
            {
                //- as long as the label of particle p is smaller
                //- than that of particle pj, collide
                //- example: p-2, pj-4,5,6
                //- so collision duplication can be avoided
                if(p.numSeq()<pj.numSeq())
                {
                    spc_.solidBinaryCollision().collide
                    (
                        p,
                        pj
                    );
                }
            }

        }
    }


    //- velocity correction
    bool noMoreCollide = false;

    for(int corStep=0;corStep<nCorrectionSteps_;corStep++)
    {
        if(!noMoreCollide)
        {
            noMoreCollide = true;
            count = 0 ;
            forAllIter(solidParticleCouplingCloud, spc_, iter)
            {
        //         solidParticleCoupling* p = &iter();
                solidParticleCoupling& p = iter();

                forAll(ngbParticles[count],j)
                {
                    //loop over neighbour particles in neighbor ngbCells
                    //- to perform binary collisions
                    solidParticleCoupling& pj = *particlePtr_()[ngbParticles[count][j]];

                    spc_.solidBinaryCollision().velocityCorrection
                    (
                        corStep,
                        p,
                        pj
                    );

                    noMoreCollide = false;
                }
            }
        }
    }
    
//     infoCounter_++;
//
//     if(infoCounter_ >= spc_.dsmcCloudReference()->nTerminalOutputs())
//     {
//             Info<< "    Inter-particle collisions       = "
//                 << collisions
//                 << endl;
//
//             infoCounter_ = 0;
//     }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
