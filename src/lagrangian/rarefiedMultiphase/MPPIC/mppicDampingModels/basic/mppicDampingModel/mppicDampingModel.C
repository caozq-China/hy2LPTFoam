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

#include "mppicDampingModel.H"
#include "mppicTimeScaleModels.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mppicDampingModel, 0);

defineRunTimeSelectionTable(mppicDampingModel, dictionary);    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class CloudType>
// Foam::mppicDampingModel<CloudType>::mppicDampingModel(CloudType& owner)
// :
//     CloudSubModelBase<CloudType>(owner),
//     timeScaleModel_(NULL)
// {}


mppicDampingModel::mppicDampingModel
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
//     const word& type
)
:
//     regionName_(dict.getOrDefault<word>("solidDampingZoneName","None")),
//     regionID_(),
    spc_(spc),
    mppicTimeScaleModel_()
{
    mppicTimeScaleModel_ = autoPtr<mppicTimeScaleModels>
    (
        mppicTimeScaleModels::New
        (
            dict.subDict("MPPICDampingModelProperties")
        )
    );
    
//     if(regionName_ != "None")
//     {
//         regionID_ = spc_.mesh().cellZones().findZoneID(regionName_);
// 
//         if(regionID_ == -1)
//         {
//             FatalErrorInFunction
//                 << "Cannot find a solid damping zone: " << regionName_ << nl << "in: "
//                 << spc_.mesh().time().system()/"mppicPropertiesDict"
//                 << exit(FatalError);
//         }
//     }
}


// template<class CloudType>
// Foam::DampingModel<CloudType>::DampingModel(const DampingModel<CloudType>& cm)
// :
//     CloudSubModelBase<CloudType>(cm),
//     timeScaleModel_(cm.timeScaleModel_)
// {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mppicDampingModel::~mppicDampingModel()
{}


// * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * * //

autoPtr<mppicDampingModel> mppicDampingModel::New
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
{
    word modelType(dict.getOrDefault<word>(typeName,"NoMppicDamping"));

    Info<< "Selecting mppicDampingModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "damping model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return
        autoPtr<mppicDampingModel>
        (
            cstrIter()(spc,dict)
        );
}

// const Foam::word& Foam::mppicDampingModel::regionName() const
// {
//     return regionName_;
// }
// 
// const Foam::labelList& Foam::mppicDampingModel::controlZone() const
// {
//     return spc_.mesh().cellZones()[regionID_];
// }


inline void Foam::mppicDampingModel::calUcorrect()
{}

} // End namespace Foam
// ************************************************************************* //
