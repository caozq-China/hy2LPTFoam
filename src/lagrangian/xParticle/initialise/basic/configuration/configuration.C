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

#include "configuration.H"
#include "xParticleCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(configuration, 0);

defineRunTimeSelectionTable(configuration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
configuration::configuration
(
    xParticleCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    initialiseDict_(dict)
    // rndGenS_(cloud.rndGenS()),
    // nSolidParticlesAdded_(0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<configuration> configuration::New
(

    xParticleCloud& cloud,
    const dictionary& dict
)
{
    const word& configurationName
    (
        dict.get<word>("type")
    );

    Info<< "Selecting configuration "
         << configurationName << endl;
         
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(configurationName);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "configuration",
            configurationName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    
    return autoPtr<configuration>
	(
		cstrIter()(cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// const word& dsmcConfiguration::name() const
// {
//     return name_;
// }

// const label& solidConfiguration::nSolidParticlesAdded() const
// {
//     return nSolidParticlesAdded_;
// }


} // End namespace Foam

// ************************************************************************* //
