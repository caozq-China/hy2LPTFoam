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

#include "absoluteMethod.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(absoluteMethod, 0);

addToRunTimeSelectionTable(mppicCorrectionLimitingMethods, absoluteMethod, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

absoluteMethod::absoluteMethod
(
    const dictionary& dict
)
:
    mppicCorrectionLimitingMethods(dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    e_(propsDict_.get<scalar>("elasticRestitutionFactor"))
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

absoluteMethod::~absoluteMethod()
{}


// * * * * * * * * * * * * * Privare Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
vector absoluteMethod::limitVelocity
(
    const vector uP,
    const vector dU,
    const vector uMean
) const
{
    const vector URelative = uP - uMean;

    return minMod
    (
        dU,
      - (1.0 + this->e_)*URelative
       *mag(uP)/max(mag(URelative), SMALL)
    );
    
}

}

// ************************************************************************* //
