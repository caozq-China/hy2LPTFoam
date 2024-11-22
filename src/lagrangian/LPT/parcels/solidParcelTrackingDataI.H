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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline Foam::solidParcel::trackingData::trackingData
(
    solidParcelCloud& cloud,
    trackPart part
)
:
    particle::TrackingData<solidParcelCloud>(cloud),
    UInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.UFilter()
        )
    ),
    rhoInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rhoFilter()
        )
    ),
    muInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.muFilter()
        )
    ),
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.TtrFilter()
        )
    ),
    CpInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.CpFilter()
        )
    ),
    kappaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.kappaFilter()
        )
    ),
    MaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.MaFilter()
        )
    ),
    gammaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.gammaFilter()
        )
    ),
    omegaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.omegaFilter()
        )
    ),
    uAverage_(),
    uSqrAverage_(),
    volumeAverage_(),
    radiusAverage_(),
    rhoAverage_(),
    frequencyAverage_(),
    massAverage_(),
    g_(cloud.g().value()),
    part_(part)
{

    if(cloud.interiorInteractionType()==2)// means MPPIC is activated
    {
        volumeAverage_ = AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":volumeAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );

        radiusAverage_ = AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":radiusAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );

        rhoAverage_ = AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":rhoAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );

        uAverage_ = AveragingMethod<vector>::New
        (
            IOobject
            (
                cloud.name() + ":uAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );

        uSqrAverage_ = AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":uSqrAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );

        frequencyAverage_ = AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":frequencyAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );

        massAverage_ = AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":massAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.subModelProperties(),
            cloud.mesh()
        );
    }
}

inline void Foam::solidParcel::trackingData::updateAverages
(
    solidParcelCloud& cloud
)
{
    // zero the sums
    volumeAverage_() = 0;
    radiusAverage_() = 0;
    rhoAverage_() = 0;
    uAverage_() = Zero;
    uSqrAverage_() = 0;
    frequencyAverage_() = 0;
    massAverage_() = 0;

    // temporary weights
    autoPtr<AveragingMethod<scalar>> weightAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                cloud.name() + ":weightAverage",
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    );
    AveragingMethod<scalar>& weightAverage = weightAveragePtr();

    // averaging sums
    forAllConstIter(solidParcelCloud, cloud, iter)
    {
        const solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        const scalar m = cloud.nParticle()*p.RWF()*p.mass();

        volumeAverage_->add(p.position(), tetIs, cloud.nParticle()*p.RWF()*p.volume());
        rhoAverage_->add(p.position(), tetIs, m*p.rho());
        uAverage_->add(p.position(), tetIs, m*p.U());
        massAverage_->add(p.position(), tetIs, m);
    }
    volumeAverage_->average();
    massAverage_->average();
    rhoAverage_->average(massAverage_);
    uAverage_->average(massAverage_);

    // squared velocity deviation
    forAllConstIter(solidParcelCloud, cloud, iter)
    {
        const solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        const vector u = uAverage_->interpolate(p.position(), tetIs);

        uSqrAverage_->add
        (
            p.position(),
            tetIs,
            cloud.nParticle()*p.RWF()*p.mass()*magSqr(p.U() - u)
        );
    }
    uSqrAverage_->average(massAverage_);

    // sauter mean radius
    radiusAverage_() = volumeAverage_();
    weightAverage = 0;
    forAllConstIter(solidParcelCloud, cloud, iter)
    {
        const solidParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        weightAverage.add
        (
            p.position(),
            tetIs,
            cloud.nParticle()*p.RWF()*pow(p.volume(), 2.0/3.0)
        );
    }
    weightAverage.average();
    radiusAverage_->average(weightAverage);

    // collision frequency
    weightAverage = 0;
    forAllConstIter(solidParcelCloud, cloud, iter)
    {
        const solidParcel& p = iter();
        tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

        const scalar a = volumeAverage_->interpolate(p.position(), tetIs);
        const scalar r = radiusAverage_->interpolate(p.position(), tetIs);
        const vector u = uAverage_->interpolate(p.position(), tetIs);

        const scalar f = 0.75*a/pow3(r)*sqr(0.5*p.d() + r)*mag(p.U() - u);

        frequencyAverage_->add(p.position(), tetIs, cloud.nParticle()*p.RWF()*f*f);

        weightAverage.add(p.position(), tetIs, cloud.nParticle()*p.RWF()*f);
        
    }
    frequencyAverage_->average(weightAverage);
}


inline const Foam::interpolation<Foam::vector>& Foam::solidParcel::trackingData::UInterp() const
{
    return UInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::rhoInterp() const
{
    return rhoInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::muInterp() const
{
    return muInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::TInterp() const
{
    return TInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::CpInterp() const
{
    return CpInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::kappaInterp() const
{
    return kappaInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::MaInterp() const
{
    return MaInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::gammaInterp() const
{
    return gammaInterp_();
}

inline const Foam::interpolation<Foam::scalar>& Foam::solidParcel::trackingData::omegaInterp() const
{
    return omegaInterp_();
}

inline const Foam::vector& Foam::solidParcel::trackingData::g() const
{
    return g_;
}

inline Foam::solidParcel::trackingData::trackPart Foam::solidParcel::trackingData::part() const
{
    return part_;
}

inline Foam::solidParcel::trackingData::trackPart& Foam::solidParcel::trackingData::part()
{
    return part_;
}


// ************************************************************************* //
