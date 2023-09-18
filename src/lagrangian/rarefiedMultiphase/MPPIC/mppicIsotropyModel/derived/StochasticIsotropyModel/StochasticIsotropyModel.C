/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "StochasticIsotropyModel.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
// class solidParticleCouplingCloud;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(StochasticIsotropyModel, 0);

addToRunTimeSelectionTable(mppicIsotropyModel, StochasticIsotropyModel, dictionary);

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class CloudType>
Foam::StochasticIsotropyModel::StochasticIsotropyModel
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    mppicIsotropyModel(spc,dict)
{}


// template<class CloudType>
// Foam::IsotropyModels::Stochastic<CloudType>::Stochastic
// (
//     const Stochastic<CloudType>& cm
// )
// :
//     IsotropyModel<CloudType>(cm)
// {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// template<class CloudType>
Foam::StochasticIsotropyModel::~StochasticIsotropyModel()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// template<class CloudType>
Foam::scalar Foam::StochasticIsotropyModel::sampleGauss()
{
    static bool isCached = true;
    static scalar xCached;

    if (isCached)
    {
        isCached = false;

        return xCached;
    }
    else
    {
//         cachedRandom& rndGen = spc_.rndGenS();

        scalar f, m, x, y;

        do
        {
            x = 2.0*spc_.rndGenS().sample01<scalar>() - 1.0;
            y = 2.0*spc_.rndGenS().sample01<scalar>() - 1.0;
            m = x*x + y*y;
        } while (m >= 1.0 || m == 0.0);

        f = sqrt(-2.0*log(m)/m);
        xCached = x*f;
        isCached = true;

        return y*f;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class CloudType>
void Foam::StochasticIsotropyModel::calculate()
{
    const fvMesh& mesh = spc_.mesh();
    const scalar deltaT(spc_.db().time().deltaTValue());
//     cachedRandom& rndGen = spc_.rndGenS();

    const scalar oneBySqrtThree = sqrt(1.0/3.0);

    const mppicAveragingMethod<scalar>& volumeAverage =
        mesh.lookupObject<mppicAveragingMethod<scalar> >
        (
            spc_.name() + ":volumeAverage"
        );
    const mppicAveragingMethod<scalar>& radiusAverage =
        mesh.lookupObject<mppicAveragingMethod<scalar> >
        (
            spc_.name() + ":radiusAverage"
        );
    const mppicAveragingMethod<vector>& uAverage =
        mesh.lookupObject<mppicAveragingMethod<vector> >
        (
            spc_.name() + ":uAverage"
        );
    const mppicAveragingMethod<scalar>& uSqrAverage =
        mesh.lookupObject<mppicAveragingMethod<scalar> >
        (
            spc_.name() + ":uSqrAverage"
        );
    const mppicAveragingMethod<scalar>& frequencyAverage =
        mesh.lookupObject<mppicAveragingMethod<scalar> >
        (
            spc_.name() + ":frequencyAverage"
        );
    const mppicAveragingMethod<scalar>& massAverage =
        mesh.lookupObject<mppicAveragingMethod<scalar> >
        (
            spc_.name() + ":massAverage"
        );

    // calculate time scales and pdf exponent
    autoPtr<mppicAveragingMethod<scalar> > exponentAveragePtr
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                spc_.name() + ":exponentAverage",
                spc_.db().time().timeName(),
                mesh
            ),
            spc_.mppicProperties(),
            mesh
        )
    );
    mppicAveragingMethod<scalar>& exponentAverage = exponentAveragePtr();
    exponentAverage =
        exp
        (
          - deltaT
           *this->mppicTimeScaleModel_->oneByTau
            (
                volumeAverage,
                radiusAverage,
                uSqrAverage,
                frequencyAverage
            )
        )();
    
    //- new version, only particles in specific zone will be corrected with return-to-isotropy model
    forAll(spc_.cellSolidCollisionOccupancyIds(), i)
    {
        const DynamicList<solidParticleCoupling*>& cellSolidParcels =  
        spc_.cellOccupancy()[spc_.cellSolidCollisionOccupancyIds()[i]];
        
        forAll(cellSolidParcels,particleID)
        {
            solidParticleCoupling& p = *cellSolidParcels[particleID];
            const tetIndices tetIs(p.currentTetIndices());

            const scalar x = exponentAverage.interpolate(p.coordinates(), tetIs);

            if (x < spc_.rndGenS().sample01<scalar>())
            {
                //- Box-Muller method
                const vector r(sampleGauss(), sampleGauss(), sampleGauss());

                const vector u = uAverage.interpolate(p.coordinates(), tetIs);
                
                //- uSqrAverage is Equation (34) in "Inclusion of collisional return-tp-isotropy in the MPPIC method" and uRms is sigma
                const scalar uRms =
                    sqrt(max(uSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));

                p.U() = u + r*uRms*oneBySqrtThree;
            }
        }
    }
    // random sampling
//     forAllIter(solidParticleCouplingCloud, spc_, iter)
//     {
//         solidParticleCoupling& p = iter();
//         
//         const tetIndices tetIs(p.currentTetIndices());
// 
//         const scalar x = exponentAverage.interpolate(p.coordinates(), tetIs);
// 
//         if (x < spc_.rndGenS().sample01<scalar>())
//         {
//             //- Box-Muller method
//             const vector r(sampleGauss(), sampleGauss(), sampleGauss());
// 
//             const vector u = uAverage.interpolate(p.coordinates(), tetIs);
//             
//             //- uSqrAverage is Equation (34) in "Inclusion of collisional return-tp-isotropy in the MPPIC method" and uRms is sigma
//             const scalar uRms =
//                 sqrt(max(uSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));
// 
//             p.U() = u + r*uRms*oneBySqrtThree;
//         }
//     }

    // correction velocity averages
    autoPtr<mppicAveragingMethod<vector> > uTildeAveragePtr
    (
        mppicAveragingMethod<vector>::New
        (
            IOobject
            (
                spc_.name() + ":uTildeAverage",
                spc_.db().time().timeName(),
                mesh
            ),
            spc_.mppicProperties(),
            mesh
        )
    );
    mppicAveragingMethod<vector>& uTildeAverage = uTildeAveragePtr();
    
//     forAllIter(solidParticleCouplingCloud, spc_, iter)
//     {
//         solidParticleCoupling& p = iter();
//         const tetIndices tetIs(p.currentTetIndices());
//         
//         scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().cellCentres()[p.cell()]);
//         
//         const scalar m = spc_.constSolidProps(p.typeID()).massSphere();
//         
//         uTildeAverage.add(p.coordinates(), tetIs, spc_.nSolidParticles()*RWF*m*p.U());
//     }
    
    forAll(spc_.cellSolidCollisionOccupancyIds(), i)
    {
        const DynamicList<solidParticleCoupling*>& cellSolidParcels =  
        spc_.cellOccupancy()[spc_.cellSolidCollisionOccupancyIds()[i]];
        
        forAll(cellSolidParcels,particleID)
        {
            solidParticleCoupling& p = *cellSolidParcels[particleID];
            const tetIndices tetIs(p.currentTetIndices());
        
            scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().cellCentres()[p.cell()]);
            
            const scalar m = spc_.constSolidProps(p.typeID()).massSphere();
            
            uTildeAverage.add(p.coordinates(), tetIs, spc_.nSolidParticles()*RWF*m*p.U());
        }
    
    }
    
    //- this is the second term on the right hand side of Equation (38) in "Inclusion of collisional return-tp-isotropy in the MPPIC method"
    uTildeAverage.average(massAverage);

    autoPtr<mppicAveragingMethod<scalar> > uTildeSqrAveragePtr
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                spc_.name() + ":uTildeSqrAverage",
                spc_.db().time().timeName(),
                mesh
            ),
            spc_.mppicProperties(),
            mesh
        )
    );
    mppicAveragingMethod<scalar>& uTildeSqrAverage = uTildeSqrAveragePtr();
    
    forAll(spc_.cellSolidCollisionOccupancyIds(), i)
    {
        const DynamicList<solidParticleCoupling*>& cellSolidParcels =  
        spc_.cellOccupancy()[spc_.cellSolidCollisionOccupancyIds()[i]];
        
        forAll(cellSolidParcels,particleID)
        {
            solidParticleCoupling& p = *cellSolidParcels[particleID];
            
            const tetIndices tetIs(p.currentTetIndices());
            const vector uTilde = uTildeAverage.interpolate(p.coordinates(), tetIs);
            
            const scalar m = spc_.constSolidProps(p.typeID()).massSphere();
            
            scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().cellCentres()[p.cell()]);
                
            uTildeSqrAverage.add
            (
                p.coordinates(),
                tetIs,
                spc_.nSolidParticles()*RWF*m*magSqr(p.U() - uTilde)
            );
        
        }
    }
//     forAllIter(solidParticleCouplingCloud, spc_, iter)
//     {
//         solidParticleCoupling& p = iter();
//         const tetIndices tetIs(p.currentTetIndices());
//         const vector uTilde = uTildeAverage.interpolate(p.coordinates(), tetIs);
//         
//         const scalar m = spc_.constSolidProps(p.typeID()).massSphere();
//         
//         scalar RWF = spc_.dsmcCloudReference()->axiRWF(spc_.mesh().cellCentres()[p.cell()]);
//             
//         uTildeSqrAverage.add
//         (
//             p.coordinates(),
//             tetIs,
//             spc_.nSolidParticles()*RWF*m*magSqr(p.U() - uTilde)
//         );
//     }

    //- this is the inverse of term on the right hand side of Equation (39) before sigmaSqr
    //- in "Inclusion of collisional return-tp-isotropy in the MPPIC method"
    uTildeSqrAverage.average(massAverage);

    //- original version 2022/11/17
    // conservation correction
//     forAllIter(solidParticleCouplingCloud, spc_, iter)
//     {
//         solidParticleCoupling& p = iter();
//         const tetIndices tetIs(p.currentTetIndices());
// 
//         const vector u = uAverage.interpolate(p.coordinates(), tetIs);
//         const scalar uRms =
//             sqrt(max(uSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));
// 
//         const vector uTilde = uTildeAverage.interpolate(p.coordinates(), tetIs);
//         const scalar uTildeRms =
//             sqrt(max(uTildeSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));
//         
//         //- this is Equation (37) in "Inclusion of collisional return-tp-isotropy in the MPPIC method"
//         p.U() = u + (p.U() - uTilde)*uRms/max(uTildeRms, SMALL);
//     }
    //- new version, only particles in specific zone will be corrected with return-to-isotropy model
    forAll(spc_.cellSolidCollisionOccupancyIds(), i)
    {
        const DynamicList<solidParticleCoupling*>& cellSolidParcels =  
        spc_.cellOccupancy()[spc_.cellSolidCollisionOccupancyIds()[i]];
        
        forAll(cellSolidParcels,particleID)
        {
            solidParticleCoupling& p = *cellSolidParcels[particleID];
            const tetIndices tetIs(p.currentTetIndices());

            const vector u = uAverage.interpolate(p.coordinates(), tetIs);
            const scalar uRms =
                sqrt(max(uSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));

            const vector uTilde = uTildeAverage.interpolate(p.coordinates(), tetIs);
            const scalar uTildeRms =
                sqrt(max(uTildeSqrAverage.interpolate(p.coordinates(), tetIs), 0.0));
            
            //- this is Equation (37) in "Inclusion of collisional return-tp-isotropy in the MPPIC method"
            p.U() = u + (p.U() - uTilde)*uRms/max(uTildeRms, SMALL);
        }
    }
    
}


bool Foam::StochasticIsotropyModel::active() const
{
    return true;
}


// ************************************************************************* //
