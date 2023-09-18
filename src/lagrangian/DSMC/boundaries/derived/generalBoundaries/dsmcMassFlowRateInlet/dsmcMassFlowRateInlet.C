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

\*---------------------------------------------------------------------------*/

#include "dsmcMassFlowRateInlet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(dsmcMassFlowRateInlet, 0);

addToRunTimeSelectionTable
(
    dsmcGeneralBoundary,
    dsmcMassFlowRateInlet,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcMassFlowRateInlet::dsmcMassFlowRateInlet
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    moleFractions_(),
    inletVelocity_(faces_.size(), Zero),
    previousInletVelocity_(faces_.size(), Zero),
    momentum_(faces_.size(), Zero),
    n_(),
    mass_(faces_.size(), 0.0),
    massFlowRate_(),
    inletTemperature_(),
    initialVelocity_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    //sorts issues with the velocity pointing out of the mesh
    inletVelocity_ = initialVelocity_;
    previousInletVelocity_ = initialVelocity_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcMassFlowRateInlet::initialConfiguration()
{}


void Foam::dsmcMassFlowRateInlet::calculateProperties()
{}


void Foam::dsmcMassFlowRateInlet::controlParcelsBeforeMove()
{
    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );

    previousInletVelocity_ = inletVelocity_;
}


void Foam::dsmcMassFlowRateInlet::controlParcelsBeforeCollisions()
{}


void Foam::dsmcMassFlowRateInlet::controlParcelsAfterCollisions()
{
    vectorField newInletVelocity(faces_.size(), Zero);

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[c]];

        for (dsmcParcel* p : parcelsInCell)
        {
            const scalar pMass =
                cloud_.nParticle()*cloud_.constProps(p->typeId()).mass();

            scalar RWF = cloud_.axiRWF(p->position());

            momentum_[c] += pMass*RWF*p->U();
            mass_[c] += pMass*RWF;
        }

        if (mass_[c] > VSMALL)
        {
            inletVelocity_[c] = momentum_[c]/mass_[c];
        }

        const vector& sF = mesh_.faceAreas()[faces_[c]];
        const scalar fA = mag(sF);

        if ((inletVelocity_[c] & -sF/fA) < 0)
        {
            inletVelocity_[c] = previousInletVelocity_[c];
        }
    }

    scalarField massFractions(typeIds_.size(), 0.0);
    scalar totalMass = 0.0;

    forAll(massFractions, iD)
    {
        const label typeId = typeIds_[iD];

        const scalar mass = cloud_.constProps(typeId).mass();

        totalMass += mass*moleFractions_[iD];
    }

    forAll(massFractions, iD)
    {
        const label typeId = typeIds_[iD];

        const scalar mass = cloud_.constProps(typeId).mass();

        massFractions[iD] = moleFractions_[iD]*(mass/totalMass);
        
        forAll(n_[iD], f)
        {
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            const scalar fA = mag(sF);
            
            n_[iD][f] = (massFractions[iD]*massFlowRate_)
               /((inletVelocity_[f] & -sF/fA)*patchSurfaceArea_*mass);
        }
    }
    
    computeParcelsToInsert
    (
        inletTemperature_,
        inletVelocity_,
        n_
    );
}


void Foam::dsmcMassFlowRateInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::dsmcMassFlowRateInlet::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    dsmcGeneralBoundary::updateProperties(dict);

    setProperties();
}


void Foam::dsmcMassFlowRateInlet::setProperties()
{
    massFlowRate_ = propsDict_.get<scalar>("massFlowRate");

    inletTemperature_ = propsDict_.get<scalar>("inletTemperature");

    initialVelocity_ = propsDict_.get<vector>("initialVelocity");

    // Read in the type ids
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // read in the mole fraction per specie

    const dictionary& moleFractionsDict = propsDict_.subDict("moleFractions");

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }

    n_.setSize(typeIds_.size());

    forAll(n_, m)
    {
        n_[m].setSize(faces_.size(), 0.0);
    }
}


// ************************************************************************* //
