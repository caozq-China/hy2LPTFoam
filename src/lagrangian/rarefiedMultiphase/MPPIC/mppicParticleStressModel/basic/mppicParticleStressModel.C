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

#include "mppicParticleStressModel.H"
// #include "solidParticleCouplingCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mppicParticleStressModel, 0);

defineRunTimeSelectionTable(mppicParticleStressModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mppicParticleStressModel::mppicParticleStressModel
(
    const dictionary& dict
)
:
    alphaPacked_(dict.get<scalar>("alphaPacked"))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<mppicParticleStressModel> mppicParticleStressModel::New
(
    const dictionary& dict
)
{
    word modelType
    (
        dict.get<word>("mppicParticleStressModelType")
    );

    Info<< "Selecting mppicParticleStressModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "particle stress model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<mppicParticleStressModel>
    (
        cstrIter()(dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mppicParticleStressModel::~mppicParticleStressModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::mppicParticleStressModel::alphaPacked() const
{
    return alphaPacked_;
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar> >
Foam::mppicParticleStressModel::tau
(
    const FieldField<Field, scalar>& alpha,
    const FieldField<Field, scalar>& rho,
    const FieldField<Field, scalar>& uRms
) const
{
    tmp<FieldField<Field, scalar> > value
    (
        new FieldField<Field, scalar>(alpha.size())
    );

    forAll(alpha, i)
    {
        value->set(i, tau(alpha[i], rho[i], uRms[i]));
    }

    return value;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
