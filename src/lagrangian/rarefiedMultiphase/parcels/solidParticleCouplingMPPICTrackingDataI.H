/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "mppicAveragingMethod.H"
#include "constants.H"


using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// #include "solidParticleCouplingCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// template<class ParcelType>
template<class CloudType>
inline Foam::solidParticleCoupling::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud
)
:
    particle::trackingData(cloud),
    volAve_
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":volumeAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    ),
    radiusAve_
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":radiusAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    ),
    rhoAve_
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":rhoAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    ),
    uAve_
    (
        mppicAveragingMethod<vector>::New
        (
            IOobject
            (
                cloud.name() + ":uAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    ),
    uSqrAve_
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":uSqrAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    ),
    frequencyAve_
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":frequencyAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    ),
    massAve_
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":massAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    )
{
}


// template<class ParcelType>
template<class CloudType>
inline void Foam::solidParticleCoupling::TrackingData<CloudType>::updateMPPICAverages
(
    CloudType& cloud
)
{
//     // zero the sums
    volAve_() = 0;
    radiusAve_() = 0;
    rhoAve_() = 0;
    uAve_() = vector::zero;
    uSqrAve_() = 0;
    frequencyAve_() = 0;
    massAve_() = 0;

    // temporary weights
    autoPtr<mppicAveragingMethod<scalar> > weightAveragePtr
    (
        mppicAveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":weightAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.mppicProperties(),
            cloud.mesh()
        )
    );
    mppicAveragingMethod<scalar>& weightAverage = weightAveragePtr();
    
    // averaging sums
//     forAllConstIter(typename CloudType, cloud, iter)
    forAllConstIter(typename CloudType, cloud, iter)
    {
        
        const solidParticleCoupling& p = iter();
        const tetIndices tetIs = p.currentTetIndices();
        
//         const scalar numSolidParticles = cloud.nSolidParticles();
        const scalar volume = pi*pow(p.D(),3.0)/6.0;
        scalar m = cloud.constSolidProps(p.typeID()).massSphere()*cloud.nSolidParticles();
        if(cloud.dsmcCloudReference()->axisymmetric())
        {
            m = m*p.RWF();
            volAve_->add(p.coordinates(), tetIs, cloud.nSolidParticles()*p.RWF()*volume);
        }
        else
        {
            volAve_->add(p.coordinates(), tetIs, cloud.nSolidParticles()*volume);
        }
        rhoAve_->add(p.coordinates(), tetIs, m*cloud.constSolidProps(p.typeID()).rho());
        uAve_->add(p.coordinates(), tetIs, m*p.U());
        massAve_->add(p.coordinates(), tetIs, m);
    }
    
    volAve_->average();
    massAve_->average();
    rhoAve_->average(*massAve_);
    uAve_->average(*massAve_);

    // squared velocity deviation
//     for (const typename CloudType::solidParticleCoupling& p : cloud)
    forAllConstIter(typename CloudType, cloud, iter)
    {
        const solidParticleCoupling& p = iter();
        const tetIndices tetIs = p.currentTetIndices();

        const vector u = uAve_->interpolate(p.coordinates(), tetIs);
        
        if(cloud.dsmcCloudReference()->axisymmetric())
        {
            uSqrAve_->add
            (
                p.coordinates(),
                tetIs,
                cloud.nSolidParticles()*p.RWF()*cloud.constSolidProps(p.typeID()).massSphere()*magSqr(p.U() - u)
            );
        }
        else
        {
            uSqrAve_->add
            (
                p.coordinates(),
                tetIs,
                cloud.nSolidParticles()*cloud.constSolidProps(p.typeID()).massSphere()*magSqr(p.U() - u)
            );
        }

        
    }
    uSqrAve_->average(*massAve_);


    // sauter mean radius
    radiusAve_() = volAve_();
    weightAverage = 0;
//     for (const typename CloudType::solidParticleCoupling& p : cloud)
    forAllConstIter(typename CloudType, cloud, iter)
    {
        
        const solidParticleCoupling& p = iter();
//         const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());
        const tetIndices tetIs = p.currentTetIndices();
        
        const scalar volume = pi*pow(p.D(),3.0)/6.0;

        if(cloud.dsmcCloudReference()->axisymmetric())
        {
            weightAverage.add
            (
                p.coordinates(),
                tetIs,
                cloud.nSolidParticles()*p.RWF()*pow(volume, 2.0/3.0)
            );
        }
        else
        {
            weightAverage.add
            (
                p.coordinates(),
                tetIs,
                cloud.nSolidParticles()*pow(volume, 2.0/3.0)
            );
        }
        
    }
    weightAverage.average();
    radiusAve_->average(weightAverage);

    // collision frequency
    weightAverage = 0;
//     for (const typename CloudType::solidParticleCoupling& p : cloud)
    forAllConstIter(typename CloudType, cloud, iter)
    {
        const solidParticleCoupling& p = iter();
//         tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());
        const tetIndices tetIs = p.currentTetIndices();

        const scalar a = volAve_->interpolate(p.coordinates(), tetIs);
        const scalar r = radiusAve_->interpolate(p.coordinates(), tetIs);
        const vector u = uAve_->interpolate(p.coordinates(), tetIs);

        const scalar f = 0.75*a/pow3(r)*sqr(0.5*p.D() + r)*mag(p.U() - u);
        
        if(cloud.dsmcCloudReference()->axisymmetric())
        {
            frequencyAve_->add(p.coordinates(), tetIs, cloud.nSolidParticles()*p.RWF()*f*f);

            weightAverage.add(p.coordinates(), tetIs, cloud.nSolidParticles()*p.RWF()*f);
        }
        else
        {
            frequencyAve_->add(p.coordinates(), tetIs, cloud.nSolidParticles()*f*f);

            weightAverage.add(p.coordinates(), tetIs, cloud.nSolidParticles()*f);
        }

        
    }
    frequencyAve_->average(weightAverage);
}



// ************************************************************************* //
