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

#include "solidFreeStreamInflowPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "solidParticleCouplingCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(solidFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable(solidGeneralBoundary, solidFreeStreamInflowPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidFreeStreamInflowPatch::solidFreeStreamInflowPatch
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidGeneralBoundary(mesh, spc, dict),
    solidPropsDict_(dict.subDict(typeName + "Properties")),
    CzRatios_()
//     phaseStates_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void solidFreeStreamInflowPatch::initialConfiguration()
{}

void solidFreeStreamInflowPatch::calculateProperties()
{}

void solidFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    computeParcelsToInsert
    (
        velocity_,
        numberDensities_
    );
    
    insertParcels
    (
        temperature_,
        velocity_,
//         phaseStates_,
        CzRatios_
    );
}

void solidFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{}

void solidFreeStreamInflowPatch::controlParcelsAfterCollisions()
{}

void solidFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void solidFreeStreamInflowPatch::updateProperties(const dictionary& Dict)
{
    //- the main properties should be updated first
    solidGeneralBoundary::updateProperties(Dict);
}



void solidFreeStreamInflowPatch::setProperties()
{

    velocity_ = solidPropsDict_.get<vector>("velocity");
    
    temperature_ = solidPropsDict_.get<scalar>("solidTemperature");

    //  read in the type ids
    typeIDs_ = spc_.getTypeIDs(solidPropsDict_);

    const dictionary& numberDensitiesSolidDict
    (
        solidPropsDict_.subDict("numberDensities")
    );
    
    //- number density initialization
    numberDensities_.clear();

    numberDensities_.setSize(typeIDs_.size(), Zero);

    forAll(numberDensities_, i)
    {
        const word& particleName = spc_.typeIdSolidList()[typeIDs_[i]];
        numberDensities_[i] = numberDensitiesSolidDict.get<scalar>(particleName);
    }

    //- cystalization ratio initialization
    CzRatios_.clear();

    CzRatios_.setSize(typeIDs_.size(), Zero);
    
    const dictionary& CzRatiosDict
    (
        solidPropsDict_.subDict("CzRatios")
    );
    
    forAll(CzRatios_, i)
    {
        const word& particleName = spc_.typeIdSolidList()[typeIDs_[i]];
        CzRatios_[i] = CzRatiosDict.get<scalar>(particleName);
        //Info<<"diametersSolid_[i] = "<<diametersSolid_[i]<<endl;
    }
    
    //- particle effective diameter initialization
//     Dps_.clear();
// 
//     Dps_.setSize(typeIDs_.size(), Zero);
    
    //- phase state initialization
//     phaseStates_.clear();
// 
//     phaseStates_.setSize(typeIDs_.size(), Zero);
//     
//     const dictionary& phaseStatesDict
//     (
//         solidPropsDict_.subDict("phaseStates")
//     );
//     
//     forAll(phaseStates_, i)
//     {
//         const word& particleName = spc_.typeIdSolidList()[typeIDs_[i]];
//         phaseStates_[i] = phaseStatesDict.get<label>(particleName);
//     }
    
    // set the accumulator
    accumulatedParcelsToInsert_.setSize(typeIDs_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(),Zero);
    }
}



} // End namespace Foam

// ************************************************************************* //
