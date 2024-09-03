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

#include "IsotropyModel.H"

#include "TimeScaleModel.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IsotropyModel, 0);

defineRunTimeSelectionTable(IsotropyModel, dictionary); 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class CloudType>
// Foam::IsotropyModel<CloudType>::IsotropyModel(CloudType& owner)
// :
//     CloudSubModelBase<CloudType>(owner),
//     timeScaleModel_(nullptr)
// {}


// template<class CloudType>
Foam::IsotropyModel::IsotropyModel
(
    const dictionary& dict,
    solidParcelCloud& cloud
    // const word& type
)
:
    cloud_(cloud),
    timeScaleModel_()
{
    timeScaleModel_ = autoPtr<TimeScaleModel>
    (
        TimeScaleModel::New
        (
            dict.subDict(TimeScaleModel::typeName)
        )
    );
}


// template<class CloudType>
// Foam::IsotropyModel<CloudType>::IsotropyModel
// (
//     const IsotropyModel<CloudType>& cm
// )
// :
//     CloudSubModelBase<CloudType>(cm),
//     timeScaleModel_(cm.timeScaleModel_)
// {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IsotropyModel::~IsotropyModel()
{}


// * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::IsotropyModel> Foam::IsotropyModel::New
(
    const dictionary& dict,
    solidParcelCloud& cloud
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting isotropy model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown isotropy model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid isotropy model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<IsotropyModel>
        (
            cstrIter()(dict, cloud)
        );
}


// ************************************************************************* //
