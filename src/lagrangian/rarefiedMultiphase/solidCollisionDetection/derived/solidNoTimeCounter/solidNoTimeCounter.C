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
    solidNoTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "solidNoTimeCounter.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidNoTimeCounter, 0);

addToRunTimeSelectionTable
(solidCollisionDetection, solidNoTimeCounter, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
solidNoTimeCounter::solidNoTimeCounter
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidCollisionDetection(mesh, spc, dict),
    infoCounter_(0)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool solidNoTimeCounter::active() const
{
    return true;
}


void solidNoTimeCounter::initialConfiguration()
{

}

void solidNoTimeCounter::collide ()
{
    if (!spc_.solidBinaryCollision().active())
    {
        return;
    }

    const scalar deltaT = spc_.dsmcCloudReference()->mesh().time().deltaTValue();

//     label collisionCandidates = 0;

    label collisions = 0;

    const List<DynamicList<solidParticleCoupling*> >& cellOccupancy = 
                                                spc_.cellOccupancy();

    const polyMesh& mesh = spc_.dsmcCloudReference()->mesh();

    forAll(cellOccupancy, cellI)
    {
        const DynamicList<solidParticleCoupling*>& cellsolidParcels(cellOccupancy[cellI]);
        const scalar cellVolume = mesh.cellVolumes()[cellI];

        label nC(cellsolidParcels.size());

        if (nC > 1)
        {
            scalar cRsigmaPPMax = spc_.cRsigmaPPMax()[cellI];
            
//             Info << "qPcRsigmaMax = " << qPcRsigmaMax << endl;
            
            scalar selectedPairs = 0.0;
               
            if(spc_.dsmcCloudReference()->axisymmetric())
            {               
                scalar RWF = 0.0;
                scalar nP = 0.0;
                
                forAll(cellsolidParcels, i)
                {
                    const solidParticleCoupling& pSolid = *cellsolidParcels[i];
                    
                    const vector position(pSolid.position());

                    scalar radius =
                        sqrt(sqr(position.y()) + sqr(position.z()));

                    RWF += 1.0 + spc_.dsmcCloudReference()->maxRWF()*(radius/spc_.dsmcCloudReference()->radialExtent());
                    
                    nP += 1.0;
                }
                
                RWF /= nP;
                
                selectedPairs = 0.5*(nC-1.0)*nC*RWF*spc_.nSolidParticles()*cRsigmaPPMax*deltaT/cellVolume;
            }
            else
            {
                
                selectedPairs = 0.5*(nC-1.0)*nC*spc_.nSolidParticles()*cRsigmaPPMax*deltaT/cellVolume;
            }

            label nCandidates(selectedPairs+spc_.rndGenS().sample01<scalar>());
            
            labelList solidParticleLabel(cellsolidParcels.size());
            
            for(label i=0; i < cellsolidParcels.size(); ++i)
            {
                solidParticleLabel[i] = i;
            }

            for (label c = 0; c < nCandidates; ++c)
            {

                // Select the first collision candidate
                label candidateP = spc_.rndGenS().position<label>(0, cellsolidParcels.size() - 1);

                // Declare the second collision candidate
                label candidateQ = -1;

                do
                {
                    candidateQ = spc_.rndGenS().position<label>(0, cellsolidParcels.size() - 1);

                } while 
                (
                    solidParticleLabel[candidateP] == solidParticleLabel[candidateQ] 
                );

                solidParticleCoupling& parcelP = *cellsolidParcels[solidParticleLabel[candidateP]];
                solidParticleCoupling& parcelQ = *cellsolidParcels[solidParticleLabel[candidateQ]];


                scalar cRsigmaPP = mag(parcelP.U()-parcelQ.U())*pi*sqr((0.5*parcelP.D()+0.5*parcelQ.D()));

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria 
                // because the number of collision candidates selected was 
                // based on this


                if (cRsigmaPP > spc_.cRsigmaPPMax()[cellI])
                {
                    spc_.cRsigmaPPMax()[cellI] = cRsigmaPP;
                }

                if ((cRsigmaPP/cRsigmaPPMax) > spc_.rndGenS().sample01<scalar>())
                {
                    //- UP != UQ

                    spc_.solidBinaryCollision().collide
                    (
                        parcelP,
                        parcelQ
                    );

                    ++collisions;
                }
                
            }
        }
    }
    

    reduce(collisions, sumOp<label>());

    spc_.cRsigmaPPMax().correctBoundaryConditions();
    
    infoCounter_++;
        
    if(infoCounter_ >= spc_.dsmcCloudReference()->nTerminalOutputs())
    {
            Info<< "    Inter-particle collisions       = "
                << collisions
                << endl;
                
            infoCounter_ = 0;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
