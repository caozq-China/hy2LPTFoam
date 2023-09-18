/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "dsmcField.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcField, 0);

defineRunTimeSelectionTable(dsmcField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dsmcField::dsmcField
(
    const Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    timeDict_(dict.subDict("timeProperties")),
    timeVel_(t, timeDict_),
    casePath_(t.path()/"fieldMeasurements"),
    timePath_()
{
    // Needed?
    mkDir(casePath_);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcField> dsmcField::New
(
    const Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("fieldModel"));

    Info<< "Selecting field: " << modelType << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "dsmcField",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<dsmcField>(cstrIter()(t, mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcField::updateTime()
{
    ++timeVel_;
}


void dsmcField::updateProperties(const dictionary& dict)
{
    timeDict_ = dict.subDict("timeProperties");

    timeDict_.readIfPresent("resetAtOutput", timeVel_.resetFieldsAtOutput());
}


const fileName& dsmcField::casePath() const
{
    return casePath_;
}


fileName& dsmcField::casePath()
{
    return casePath_;
}


const fileName& dsmcField::timePath() const
{
    return timePath_;
}


fileName& dsmcField::timePath()
{
    return timePath_;
}

} // End namespace Foam

// ************************************************************************* //
