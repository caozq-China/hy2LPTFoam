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

#include "RelaxationDamping.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// class solidParticleCouplingCloud;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RelaxationDamping, 0);

addToRunTimeSelectionTable(mppicDampingModel, RelaxationDamping, dictionary);

}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RelaxationDamping::RelaxationDamping
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    mppicDampingModel(spc, dict),
    uAverage_(nullptr),
    oneByTimeScaleAverage_(nullptr)
{
    
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::RelaxationDamping::~RelaxationDamping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class CloudType>
void Foam::RelaxationDamping::cacheFields(const bool store)
{
    if (store)
    {
        const fvMesh& mesh = spc_.mesh();
        const word& cloudName = spc_.name();

        const mppicAveragingMethod<scalar>& volumeAverage =
            mesh.lookupObject<mppicAveragingMethod<scalar> >
            (
                cloudName + ":volumeAverage"
            );
        const mppicAveragingMethod<scalar>& radiusAverage =
            mesh.lookupObject<mppicAveragingMethod<scalar> >
            (
                cloudName + ":radiusAverage"
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
        const mppicAveragingMethod<scalar>& frequencyAverage =
            mesh.lookupObject<mppicAveragingMethod<scalar> >
            (
                cloudName + ":frequencyAverage"
            );

        uAverage_ = &uAverage;

        oneByTimeScaleAverage_.reset
        (
            mppicAveragingMethod<scalar>::New
            (
                IOobject
                (
                    cloudName + ":oneByTimeScaleAverage",
                    spc_.db().time().timeName(),
                    mesh
                ),
                spc_.mppicProperties(),
                mesh
            ).ptr()
        );

        oneByTimeScaleAverage_() =
        (
            this->mppicTimeScaleModel_->oneByTau
            (
                volumeAverage,
                radiusAverage,
                uSqrAverage,
                frequencyAverage
            )
        )();
    }
    else
    {
        uAverage_ = nullptr;
        oneByTimeScaleAverage_.clear();
    }
}


// template<class CloudType>
Foam::vector Foam::RelaxationDamping::velocityCorrection
(
    solidParticleCoupling& p,
    const scalar deltaT
) const
{
    const tetIndices tetIs(p.currentTetIndices());

    const scalar x =
        deltaT*oneByTimeScaleAverage_->interpolate(p.coordinates(), tetIs);

    const vector u = uAverage_->interpolate(p.coordinates(), tetIs);

    return (u - p.U())*x/(x + 2.0);
}

void Foam::RelaxationDamping::calUcorrect()
{
//     Info<<"Enter Damping Zone Correction"<<endl;
    
//     const label zoneId = spc_.mesh().cellZones().findZoneID(regionName_);
//     const labelList& zone = spc_.mesh().cellZones()[zoneId];
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
//                 pSolid.UCorrect() = spc_.mppicDampingModels().velocityCorrection(pSolid,deltaT);
//             }
//         }
//     }
//     else
//     {
//         forAllIter(solidParticleCouplingCloud, spc, iter)
//         {
//             solidParticleCoupling* pSolid = &iter();
//             
//             pSolid->UCorrect() = spc_.mppicDampingModels().velocityCorrection(iter(),deltaT);
//         }
//     }
    forAll(spc_.cellSolidPackingOccupancyIds(),i)
    {
        const DynamicList<solidParticleCoupling*>& cellSolidParcels =  spc_.cellOccupancy()[spc_.cellSolidPackingOccupancyIds()[i]];
        
        forAll(cellSolidParcels,particleID)
        {
            solidParticleCoupling& pSolid = *cellSolidParcels[particleID];
            pSolid.UCorrect() = spc_.mppicDampingModels().velocityCorrection(pSolid,spc_.mesh().time().deltaTValue());
        }
    }
    
}


bool Foam::RelaxationDamping::active() const
{
    return true;
}


// ************************************************************************* //
