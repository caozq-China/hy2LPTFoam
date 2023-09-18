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

#include "ExplicitPacking.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ExplicitPacking, 0);

addToRunTimeSelectionTable(mppicPackingModel, ExplicitPacking, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ExplicitPacking::ExplicitPacking
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    mppicPackingModel(spc,dict),
    velocityAve_(),
    volumeAverage_(),
    stressAverage_(nullptr),
    mppicCorrectionLimitingMethods_()
{
    
    mppicCorrectionLimitingMethods_ = autoPtr<mppicCorrectionLimitingMethods>
    (
        mppicCorrectionLimitingMethods::New
        (
            dict
        )
    );
    
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ExplicitPacking::~ExplicitPacking()
{}


// * * * * * * * * * * * * * Privare Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::ExplicitPacking::cacheFields(const bool store)
{
    mppicPackingModel::cacheFields(store);
    
    if(store)
    {
        const fvMesh& mesh = spc_.mesh();
        const word& cloudName = spc_.name();

        const mppicAveragingMethod<scalar>& volumeAverage =
            mesh.lookupObject<mppicAveragingMethod<scalar> >
            (
                cloudName + ":volumeAverage"
            );
        const mppicAveragingMethod<scalar>& rhoAverage =
            mesh.lookupObject<mppicAveragingMethod<scalar> >
            (
                cloudName + ":rhoAverage"
            );
        const mppicAveragingMethod<vector>& uAverage =
            mesh.lookupObject<mppicAveragingMethod<vector> >
            (
                cloudName + ":uAverage"
            );
        const mppicAveragingMethod<scalar>& uSqrAverage =
            mesh.lookupObject<mppicAveragingMethod<scalar> >
            (
                cloudName + ":uSqrAverage"
            );

        volumeAverage_ = &volumeAverage;
        velocityAve_ = &uAverage;

        stressAverage_.reset
        (
            mppicAveragingMethod<scalar>::New
            (
                IOobject
                (
                    cloudName + ":stressAverage",
                    spc_.db().time().timeName(),
                    mesh
                ),
                spc_.mppicProperties(),
                mesh
            ).ptr()
        );
        stressAverage_() =
            this->mppicParticleStressModel_->tau
            (
                *volumeAverage_,
                rhoAverage,
                uSqrAverage
            )();
    }
    else
    {
        volumeAverage_ = nullptr;
        velocityAve_ = nullptr;
        stressAverage_.clear();
    }
}

Foam::vector Foam::ExplicitPacking::velocityUpdateAndCorrection
(
    solidParticleCoupling& p,
    const scalar deltaT
) const
{
//     const fvMesh& mesh = spc_.mesh();

    const tetIndices tetIs(p.currentTetIndices());
    // interpolated quantities
    const scalar alpha =
        this->volumeAverage_->interpolate(p.coordinates(), tetIs);
    const vector alphaGrad =
        this->volumeAverage_->interpolateGrad(p.coordinates(), tetIs);
    const vector uMean =
        this->velocityAve_->interpolate(p.coordinates(), tetIs);

    // stress gradient
    const vector tauGrad =
        stressAverage_->interpolateGrad(p.coordinates(), tetIs);

    // parcel relative velocity
    const vector uRelative = p.U() - uMean;

    // correction velocity
    vector dU = Zero;

    // correction velocity
    if ((uRelative & alphaGrad) > 0)
    {
        dU = - deltaT*tauGrad/(spc_.constSolidProps(p.typeID()).rho()*(alpha + SMALL));
    }

    // apply the velocity limiters
    return
        mppicCorrectionLimitingMethods_->limitVelocity
        (
            p.U(),
            dU,
            uMean
        );
//     return vector::zero;
}

void Foam::ExplicitPacking::calUcorrect()
{

//     if(regionName() != "None")
//     {
//         forAll(controlZone(), c)
//         {
//     //         const auto& cellOccupancy = spc_.cellOccupancy();
//             const label cellI = controlZone()[c];
//             const DynamicList<solidParticleCoupling*>& cellSolidParcels =  spc_.cellOccupancy()[cellI];
//             
//             forAll(cellSolidParcels,particleID)
//             {
//                 solidParticleCoupling& pSolid = *cellSolidParcels[particleID];
//                 pSolid.UCorrect() = spc_.mppicPackingModels().velocityUpdateAndCorrection(pSolid,deltaT);
//             }
//         }
//     }
//     else
//     {
//         forAllIter(solidParticleCouplingCloud, spc, iter)
//         {
//             solidParticleCoupling* pSolid = &iter();
//             
//             pSolid->UCorrect() = spc_.mppicPackingModels().velocityUpdateAndCorrection(iter(),deltaT);
//         }
//     }
    forAll(spc_.cellSolidPackingOccupancyIds(),i)
    {
        const DynamicList<solidParticleCoupling*>& cellSolidParcels =  spc_.cellOccupancy()[spc_.cellSolidPackingOccupancyIds()[i]];
        
        forAll(cellSolidParcels,particleID)
        {
            solidParticleCoupling& pSolid = *cellSolidParcels[particleID];
            pSolid.UCorrect() = spc_.mppicPackingModels().velocityUpdateAndCorrection(pSolid,spc_.mesh().time().deltaTValue());
        }
    }
    
}


bool Foam::ExplicitPacking::active() const
{
    return true;
}

}//End namespace Foam

// ************************************************************************* //
