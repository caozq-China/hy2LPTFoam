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
    collisionPartnerSelection

Description

\*----------------------------------------------------------------------------*/

#include "mppicPackingModel.H"
#include "mppicAveragingMethod.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mppicPackingModel, 0);

defineRunTimeSelectionTable(mppicPackingModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mppicPackingModel::mppicPackingModel
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
):
//     regionName_(dict.getOrDefault<word>("solidPackingZoneName","None")),
//     regionID_(),
    spc_(spc),
    mppicParticleStressModel_()
{
    mppicParticleStressModel_ = 
        mppicParticleStressModel::New
        (
            dict
        );
//     if(regionName_ != "None")
//     {
//         regionID_ = spc_.mesh().cellZones().findZoneID(regionName_);
// 
//         if(regionID_ == -1)
//         {
//             FatalErrorInFunction
//                 << "Cannot find a solid packing zone: " << regionName_ << nl << "in: "
//                 << spc_.mesh().time().system()/"mppicPropertiesDict"
//                 << exit(FatalError);
//         }
//     }
}


// mppicPackingModel::mppicPackingModel(const mppicPackingModel& PM)
// :
//     (PM),
//     mppicParticleStressModel_(PM.mppicParticleStressModel_)
// {}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<mppicPackingModel> mppicPackingModel::New
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
{
    word modelType
    (
        dict.getOrDefault<word>(typeName, "NoPacking")
    );

    Info<< "Selecting mppicPackingModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "packing model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<mppicPackingModel>
    (
        cstrIter()(spc, dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mppicPackingModel::~mppicPackingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// const Foam::word& Foam::mppicPackingModel::regionName() const
// {
//     return regionName_;
// }

// const Foam::labelList& Foam::mppicPackingModel::controlZone() const
// {
//     return spc_.mesh().cellZones()[regionID_];
// }

inline void Foam::mppicPackingModel::cacheFields(const bool store)
{}

inline void Foam::mppicPackingModel::calUcorrect(){}

inline const Foam::mppicParticleStressModel& Foam::mppicPackingModel::mppicParticleStress() const
{
    return mppicParticleStressModel_;
}

inline Foam::mppicParticleStressModel& Foam::mppicPackingModel::mppicParticleStress()
{
    return mppicParticleStressModel_();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
