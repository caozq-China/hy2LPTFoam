/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "Stochastic.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Stochastic, 0);

addToRunTimeSelectionTable(IsotropyModel, Stochastic, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Stochastic::Stochastic
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
:
    IsotropyModel(dict, cloud)
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

Foam::Stochastic::~Stochastic()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// template<class CloudType>
Foam::scalar Foam::Stochastic::sampleGauss()
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
        // Random& rndGen = this->owner().rndGen();

        scalar f, m, x, y;

        do
        {
            x = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            y = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
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
void Foam::Stochastic::calculate()
{
    const fvMesh& mesh = cloud_.mesh();
    const scalar deltaT(cloud_.db().time().deltaTValue());
    // Random& rndGen = cloud_.rndGen();

    const scalar oneBySqrtThree = sqrt(1.0/3.0);

    const AveragingMethod<scalar>& volumeAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            cloud_.name() + ":volumeAverage"
        );
    const AveragingMethod<scalar>& radiusAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            cloud_.name() + ":radiusAverage"
        );
    const AveragingMethod<vector>& uAverage =
        mesh.lookupObject<AveragingMethod<vector>>
        (
            cloud_.name() + ":uAverage"
        );
    const AveragingMethod<scalar>& uSqrAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            cloud_.name() + ":uSqrAverage"
        );
    const AveragingMethod<scalar>& frequencyAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            cloud_.name() + ":frequencyAverage"
        );
    const AveragingMethod<scalar>& massAverage =
        mesh.lookupObject<AveragingMethod<scalar>>
        (
            cloud_.name() + ":massAverage"
        );

    // calculate time scales and pdf exponent
    autoPtr<AveragingMethod<scalar>> exponentAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud_.name() + ":exponentAverage",
                cloud_.db().time().timeName(),
                mesh
            ),
            cloud_.subModelProperties(),
            mesh
        )
    );
    AveragingMethod<scalar>& exponentAverage = exponentAveragePtr();
    exponentAverage =
        exp
        (
          - deltaT
           *this->timeScaleModel_->oneByTau
            (
                volumeAverage,
                radiusAverage,
                uSqrAverage,
                frequencyAverage
            )
        )();

    // random sampling
    forAllIter(solidParcelCloud, cloud_, iter)
    {
        solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh);

        const scalar x = exponentAverage.interpolate(p.position(), tetIs);

        if (x < cloud_.rndGen().sample01<scalar>())
        {
            const vector r(sampleGauss(), sampleGauss(), sampleGauss());

            const vector u = uAverage.interpolate(p.position(), tetIs);
            const scalar uRms =
                sqrt(max(uSqrAverage.interpolate(p.position(), tetIs), 0.0));

            p.U() = u + r*uRms*oneBySqrtThree;
        }
    }

    // correction velocity averages
    autoPtr<AveragingMethod<vector>> uTildeAveragePtr
    (
        AveragingMethod<vector>::New
        (
            IOobject
            (
                cloud_.name() + ":uTildeAverage",
                cloud_.db().time().timeName(),
                mesh
            ),
            cloud_.subModelProperties(),
            mesh
        )
    );
    AveragingMethod<vector>& uTildeAverage = uTildeAveragePtr();
    forAllIter(solidParcelCloud, cloud_, iter)
    {
        solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh);
        // uTildeAverage.add(p.position(), tetIs, cloud_.nParticle()*p.mass()*p.U());
        uTildeAverage.add(p.position(), tetIs, cloud_.nParticle()*p.RWF()*p.mass()*p.U());
    }
    uTildeAverage.average(massAverage);

    autoPtr<AveragingMethod<scalar>> uTildeSqrAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud_.name() + ":uTildeSqrAverage",
                cloud_.db().time().timeName(),
                mesh
            ),
            cloud_.subModelProperties(),
            mesh
        )
    );
    AveragingMethod<scalar>& uTildeSqrAverage = uTildeSqrAveragePtr();
    forAllIter(solidParcelCloud, cloud_, iter)
    {
        solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh);
        const vector uTilde = uTildeAverage.interpolate(p.position(), tetIs);
        uTildeSqrAverage.add
        (
            p.position(),
            tetIs,
            cloud_.nParticle()*p.RWF()*p.mass()*magSqr(p.U() - uTilde)
            // cloud_.nParticle()*p.mass()*magSqr(p.U() - uTilde)
        );
    }
    uTildeSqrAverage.average(massAverage);

    // conservation correction
    forAllIter(solidParcelCloud, cloud_, iter)
    {
        solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh);

        const vector u = uAverage.interpolate(p.position(), tetIs);
        const scalar uRms =
            sqrt(max(uSqrAverage.interpolate(p.position(), tetIs), 0.0));

        const vector uTilde = uTildeAverage.interpolate(p.position(), tetIs);
        const scalar uTildeRms =
            sqrt(max(uTildeSqrAverage.interpolate(p.position(), tetIs), 0.0));

        p.U() = u + (p.U() - uTilde)*uRms/max(uTildeRms, SMALL);
    }
}

bool Foam::Stochastic::active() const
{
    return true;
}

// ************************************************************************* //
