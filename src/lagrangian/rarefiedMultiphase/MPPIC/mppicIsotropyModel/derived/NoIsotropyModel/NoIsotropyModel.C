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

#include "NoIsotropyModel.H"
#include "addToRunTimeSelectionTable.H"



namespace Foam
{
// class solidParticleCouplingCloud;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NoIsotropyModel, 0);

addToRunTimeSelectionTable(mppicIsotropyModel, NoIsotropyModel, dictionary);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class CloudType>
Foam::NoIsotropyModel::NoIsotropyModel
(
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    mppicIsotropyModel(spc,dict)
{}


// template<class CloudType>
// Foam::IsotropyModels::NoIsotropy<CloudType>::NoIsotropy
// (
//     const NoIsotropy<CloudType>& cm
// )
// :
//     IsotropyModel<CloudType>(cm)
// {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// template<class CloudType>
Foam::NoIsotropyModel::~NoIsotropyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class CloudType>
void Foam::NoIsotropyModel::calculate()
{
    // do nothing
}


// template<class CloudType>
bool Foam::NoIsotropyModel::active() const
{
    return false;
}


// ************************************************************************* //
