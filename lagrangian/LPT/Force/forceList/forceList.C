/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "forceList.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::forceList::forceList
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool readFields
)
:
    mesh_(mesh),
    dict_(dict),
    calcCoupled_(true),
    calcNonCoupled_(true)
{
    if (readFields)
    {
        wordList modelNames(dict.toc());

        Info<< "Constructing particle forces" << endl;

        if (modelNames.size() > 0)
        {
            this->setSize(modelNames.size());

            label i = 0;
            forAllConstIter(IDLList<entry>, dict, iter)
            {
                const word& model = iter().keyword();

                if (iter().isDict())
                {
                    // Info<<"1"<<endl;
                    this->set
                    (
                        i++,
                        Force::New
                        (
                            cloud,
                            mesh,
                            iter().dict(),
                            model
                        )
                    );
                }
                else
                {
                    // Info<<"2"<<endl;
                    this->set
                    (
                        i++,
                        Force::New
                        (
                            cloud,
                            mesh,
                            dict,
                            model
                        )
                    );
                }
            }
        }
        else
        {
            Info<< "    none" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::forceList::~forceList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forceList::cacheFields(const bool store)
{
    forAll(*this, i)
    {
        this->operator[](i).cacheFields(store);
    }
}

Foam::forceSuSp Foam::forceList::calcCoupled
(
    const Foam::solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    if(calcCoupled_)
    {
        forAll(*this, i)
        {
            value += this->operator[](i).calcCoupled(p, dt, mass, Re, muc);
        }
    }

    return value;
}

Foam::forceSuSp Foam::forceList::calcNonCoupled
(
    const Foam::solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    if(calcNonCoupled_)
    {
        forAll(*this, i)
        {
            value += this->operator[](i).calcNonCoupled(p, dt, mass, Re, muc);
        }
    }
    return value;
}

Foam::scalar Foam::forceList::massEff
(
    const solidParcel& p,
    const scalar mass
) const
{
    scalar massEff = mass;
    forAll(*this, i)
    {
        massEff += this->operator[](i).massAdd(p, mass);
    }

    return massEff;
}


// ************************************************************************* //
