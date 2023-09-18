/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "solidField.H"
// #include "graph.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidField, 0);

defineRunTimeSelectionTable(solidField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidField::solidField
(
    const Time& t,
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
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
    mkDir(casePath_);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidField> solidField::New
(
    const Time& t,
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
{
    const word modelType
    (
        dict.get<word>("fieldModel")
    );

    Info<< "Selecting field: "
         << modelType << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);
    
    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidField",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<solidField>
    (
        cstrIter()(t, mesh, cloud, dict)
    );
}

void solidField::updateTime()
{
    ++timeVel_;
}


void solidField::updateProperties
(
    const dictionary& dict
)
{
    timeDict_ = dict.subDict("timeProperties");
    
    timeDict_.readIfPresent("resetAtOutput", timeVel_.resetFieldsAtOutput());

}

const fileName& solidField::casePath() const
{
    return casePath_;
}

fileName& solidField::casePath()
{
    return casePath_;
}


const fileName& solidField::timePath() const
{
    return timePath_;
}

fileName& solidField::timePath()
{
    return timePath_;
}



} // End namespace Foam

// ************************************************************************* //
