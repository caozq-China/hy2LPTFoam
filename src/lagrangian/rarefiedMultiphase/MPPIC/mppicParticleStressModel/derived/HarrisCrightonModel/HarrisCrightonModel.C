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

#include "HarrisCrightonModel.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(HarrisCrightonModel, 0);

addToRunTimeSelectionTable(mppicParticleStressModel, HarrisCrightonModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HarrisCrightonModel::HarrisCrightonModel
(
    const dictionary& dict
)
:
    mppicParticleStressModel(dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    pressureCoeff_(propsDict_.get<scalar>("pressureCoeff")),
    beta_(propsDict_.get<scalar>("beta")),
    eps_(propsDict_.get<scalar>("eps"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
HarrisCrightonModel::~HarrisCrightonModel()
{}


// * * * * * * * * * * * * * Privare Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar> > Foam::HarrisCrightonModel::denominator
(
    const Field<scalar>& alpha
) const
{
    return
        max
        (
            alphaPacked_ - alpha,
            max(eps_*(1.0 - alpha), SMALL)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar> > Foam::HarrisCrightonModel::tau
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    return
    (
        pressureCoeff_
      * pow(alpha, beta_)
      / denominator(alpha)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::HarrisCrightonModel::dTaudTheta
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    const Field<scalar> d(denominator(alpha));
    
    return
    (
        pressureCoeff_
      * pow(alpha, beta_)
      / d
      * (beta_/alpha + 1.0/d)
    );
}


}

// ************************************************************************* //
