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

#include "solidConfiguration.H"
// #include "IFstream.H"
// #include "graph.H"
#include "solidParticleCouplingCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidConfiguration, 0);

defineRunTimeSelectionTable(solidConfiguration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidConfiguration::solidConfiguration
(
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    solidInitialiseDict_(dict),
    rndGenS_(cloud.rndGenS()),
    nSolidParticlesAdded_(0)
{

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidConfiguration> solidConfiguration::New
(

    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
{
    const word& solidConfigurationName
    (
        dict.get<word>("type")
    );

    Info<< "Selecting solidConfiguration "
         << solidConfigurationName << endl;
         
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(solidConfigurationName);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidConfiguration",
            solidConfigurationName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    
    return autoPtr<solidConfiguration>
	(
		cstrIter()(cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// solidConfiguration::~solidConfiguration()
// {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// vector dsmcConfiguration::equipartitionLinearVelocity
// (
//     scalar temperature,
//     scalar mass
// )
// {
//     return sqrt(molCloud_.redUnits().kB()*temperature/mass)*vector
//     (
//         rndGen_.GaussNormal(),
//         rndGen_.GaussNormal(),
//         rndGen_.GaussNormal()
//     );
// }

// const word& dsmcConfiguration::name() const
// {
//     return name_;
// }

const label& solidConfiguration::nSolidParticlesAdded() const
{
    return nSolidParticlesAdded_;
}


} // End namespace Foam

// ************************************************************************* //
