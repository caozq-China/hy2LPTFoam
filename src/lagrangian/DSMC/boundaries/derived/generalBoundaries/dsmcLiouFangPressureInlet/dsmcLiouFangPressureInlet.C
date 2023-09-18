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

#include "dsmcLiouFangPressureInlet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(dsmcLiouFangPressureInlet, 0);

addToRunTimeSelectionTable
(
    dsmcGeneralBoundary,
    dsmcLiouFangPressureInlet,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dsmcLiouFangPressureInlet::dsmcLiouFangPressureInlet
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
    inletPressure_(),
    inletTemperature_(),
    theta_(),
    n_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    n_ = inletPressure_ / (physicoChemical::k.value()*inletTemperature_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcLiouFangPressureInlet::initialConfiguration()
{}


void dsmcLiouFangPressureInlet::calculateProperties()
{}


void dsmcLiouFangPressureInlet::controlParcelsBeforeMove()
{
    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );

    previousInletVelocity_ = inletVelocity_;
}


void dsmcLiouFangPressureInlet::controlParcelsBeforeCollisions()
{}


void dsmcLiouFangPressureInlet::controlParcelsAfterCollisions()
{
    vectorField momentum(faces_.size(), Zero);
    vectorField newInletVelocity(faces_.size(), Zero);
    scalarField mass(faces_.size(), scalar(0));


    const List<DynamicList<dsmcParcel*>>& cellOccupancy =
        cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[c]];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];

            const scalar m =
                cloud_.nParticle()*cloud_.constProps(p->typeId()).mass();

            scalar RWF = cloud_.axiRWF(p->position());

            momentum[c] += RWF*m*p->U();
            mass[c] += RWF*m;
        }

        newInletVelocity[c] = momentum[c]/mass[c];

        inletVelocity_[c] =
            theta_*newInletVelocity[c]
          + (1.0 - theta_)*previousInletVelocity_[c];
    }

    // Compute number of parcels to insert
    computeParcelsToInsert
    (
        inletTemperature_,
        inletVelocity_,
        n_,
        moleFractions_
    );
}


void dsmcLiouFangPressureInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcLiouFangPressureInlet::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    dsmcGeneralBoundary::updateProperties(dict);

    setProperties();
}


void dsmcLiouFangPressureInlet::setProperties()
{
    inletPressure_ = propsDict_.get<scalar>("inletPressure");

    inletTemperature_ = propsDict_.get<scalar>("inletTemperature");

    theta_ = propsDict_.get<scalar>("theta");

    if (0.0 > theta_ || theta_ > 1.0)
    {
        FatalErrorInFunction
            << "Theta must be a value between 0 and 1 " << nl << "in: "
            << mesh_.time().system()/dsmcBoundaries::dictName
            << exit(FatalError);
    }

    // Read in the type ids
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Read in the mole fraction per specie

    const dictionary& moleFractionsDict
    (
        propsDict_.subDict("moleFractions")
    );

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    // Set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }
}

} // End namespace Foam

// ************************************************************************* //
