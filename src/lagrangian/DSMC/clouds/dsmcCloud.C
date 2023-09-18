/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "dsmcCloud.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcCloud::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] = dsmcParcel::constantProperties(molDict);
    }
}


void Foam::dsmcCloud::removeElectrons()
{
    const auto& cellOccupancy = this->cellOccupancy();

    forAll(cellOccupancy, celli)
    {
        const auto& occupancy = cellOccupancy[celli];

        for (dsmcParcel* p : occupancy)
        {
            const auto& constProp = constProps(p->typeId());

            scalar RWF = axiRWF(mesh_.cellCentres()[celli]);

            scalar mass = constProps(p->typeId()).mass();

            momentumMean_[celli] += mass*RWF*p->U();
            rhoMMean_[celli] += mass*RWF;

            if (constProp.charge() == -1)
            {
                rhoNMeanElectron_[celli] += 1.0*RWF;
                rhoMMeanElectron_[celli] += mass*RWF;
                momentumMeanElectron_[celli] += mass*RWF*p->U();
                linearKEMeanElectron_[celli] += mass*RWF*(p->U() & p->U());

                // found an electron
                deleteParticle(*p);
            }
        }
    }
}


void Foam::dsmcCloud::addElectrons()
{
    label electronTypeId = -1;

    // find electron typeId
    forAll(constProps_, cP)
    {
        const label electronCharge = constProps_[cP].charge();

        if (electronCharge == -1)
        {
            electronTypeId = cP;
            break;
        }
    }

    auto& cellOccupancy = this->cellOccupancy();

    forAll(cellOccupancy, celli)
    {
        if (rhoMMeanElectron_[celli] > VSMALL)
        {
            scalar V = mesh_.cellVolumes()[celli];

            scalar rhoMMeanElectron = rhoMMeanElectron_[celli]*nParticle_/V;
            scalar rhoNMeanElectron = rhoNMeanElectron_[celli]*nParticle_/V;
            vector UElectron =
                momentumMeanElectron_[celli]/(rhoMMeanElectron*V);
            scalar linearKEMeanElectron =
                (0.5*linearKEMeanElectron_[celli]*nParticle_)/V;

            electronTemperature_[celli] =
                2.0
               /(3.0*physicoChemical::k.value()*rhoNMeanElectron)
               *(
                    linearKEMeanElectron
                  - 0.5*rhoMMeanElectron*(UElectron & UElectron)
                );
        }

        for (dsmcParcel* p : cellOccupancy[celli])
        {
            const auto& constProp = constProps(p->typeId());

            if (constProp.charge() == 1)
            {
                // found an ion, add an electron here

                // electron temperature will be zero if there have been no
                // electrons in the cell during the simulation


                const vector position(p->position());

                if (electronTemperature_[celli] < VSMALL)
                {
                    electronTemperature_[celli] = 6000.0;
                }
                if (electronTemperature_[celli] > 8.0e4)
                {
                    electronTemperature_[celli] = 30000.0;
                }

                vector electronVelocity =
                    equipartitionLinearVelocity
                    (
                        electronTemperature_[celli],
                        constProps_[electronTypeId].mass()
                    );

                if (rhoMMean_[celli] > VSMALL)
                {
                    cellVelocity_[celli] =
                        momentumMean_[celli]/rhoMMean_[celli];
                }

                labelList vibLevel(0, 0);

                electronVelocity += cellVelocity_[celli];

                scalar RWF = p->RWF();

                addNewParcel
                (
                    position,
                    electronVelocity,
                    RWF,
                    0.0,
                    0,
                    celli,
                    electronTypeId,
                    0,
                    vibLevel
                );
            }
        }
    }
}


Foam::label Foam::dsmcCloud::pickFromCandidateList
(
    DynamicList<label>& candidatesInCell
)
{
    label entry = -1;
    const label size = candidatesInCell.size();

    if (size > 0)
    {
        // choose a random number between 0 and the size of the candidateList
        // size
        label randomIndex = rndGen_.position<label>(0, size-1);
        entry = candidatesInCell[randomIndex];

        // build a new list without the chosen entry
        DynamicList<label> newCandidates(0);

        forAll(candidatesInCell, i)
        {
            if (i != randomIndex)
            {
                newCandidates.append(candidatesInCell[i]);
            }
        }

        // transfer the new list
        candidatesInCell.transfer(newCandidates);
        candidatesInCell.shrink();
    }

    return entry;
}


void Foam::dsmcCloud::updateCandidateSubList
(
    const label candidate,
    DynamicList<label>& candidatesInSubCell
)
{
    const label newIndex = candidatesInSubCell.find(candidate);

    DynamicList<label> newCandidates(0);

    forAll(candidatesInSubCell, i)
    {
        if (i != newIndex)
        {
            newCandidates.append(candidatesInSubCell[i]);
        }
    }

    // transfer the new list
    candidatesInSubCell.transfer(newCandidates);
    candidatesInSubCell.shrink();
}


Foam::label Foam::dsmcCloud::pickFromCandidateSubList
(
    DynamicList<label>& candidatesInCell,
    DynamicList<label>& candidatesInSubCell
)
{
    label entry = -1;
    const label subCellSize = candidatesInSubCell.size();

    if (subCellSize > 0)
    {
        label randomIndex = rndGen_.position<label>(0, subCellSize-1);
        entry = candidatesInSubCell[randomIndex];

        DynamicList<label> newSubCellList(0);

        forAll(candidatesInSubCell, i)
        {
            if (i != randomIndex)
            {
                newSubCellList.append(candidatesInSubCell[i]);
            }
        }

        candidatesInSubCell.transfer(newSubCellList);
        candidatesInSubCell.shrink();

        label newIndex = candidatesInCell.find(entry);

        DynamicList<label> newList(0);

        forAll(candidatesInCell, i)
        {
            if (i != newIndex)
            {
                newList.append(candidatesInCell[i]);
            }
        }

        candidatesInCell.transfer(newList);
        candidatesInCell.shrink();
    }

    return entry;
}


void Foam::dsmcCloud::collisions()
{
    collisionPartnerSelectionModel_->collide();
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::dsmcCloud::addNewParcel
(
    const vector& position,
    const vector& U,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label typeId,
    const label newParcel,
    const labelList vibLevel
)
{
    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        RWF,
        ERot,
        ELevel,
        cellI,
        typeId,
        newParcel,
        vibLevel
    );

    addParticle(pPtr);
}


Foam::scalar Foam::dsmcCloud::energyRatio
(
    const scalar ChiA,
    const scalar ChiB
)
{
    const scalar ChiAMinusOne = ChiA - 1;
    const scalar ChiBMinusOne = ChiB - 1;

    if (ChiAMinusOne<SMALL && ChiBMinusOne<SMALL)
    {
        return rndGen_.sample01<scalar>();
    }

    scalar energyRatio;

    scalar P;

    const scalar eps = rndGen_.sample01<scalar>();

    do
    {
        P = 0;

        energyRatio = rndGen_.sample01<scalar>();

        if (ChiAMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio), ChiBMinusOne);
        }
        else if (ChiBMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio), ChiAMinusOne);
        }
        else
        {
            P =
                pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                    ChiAMinusOne
                )
               *pow
                (
                    (ChiAMinusOne + ChiBMinusOne)
                   *(1 - energyRatio)
                   /ChiBMinusOne,
                    ChiBMinusOne
                );
        }
    } while (P < eps);

    return energyRatio;
}


Foam::scalar Foam::dsmcCloud::PSIm
(
    const scalar DOFm,
    const scalar DOFtot
)
{
    if (DOFm == DOFtot)
    {
        return 1.0;
    }

    if (DOFm == 2.0 && DOFtot == 4.0)
    {
        return rndGen_.sample01<scalar>();
    }

    if (DOFtot < 4.0)
    {
        return (DOFm/DOFtot);
    }

    scalar rPSIm = 0.0;
    scalar prob = 0.0;

    const scalar h1 = 0.5*DOFtot - 2.0;
    const scalar h2 = 0.5*DOFm - 1.0 + 1.0e-5;
    const scalar h3 = 0.5*(DOFtot - DOFm)-1.0 + 1.0e-5;

    const scalar eps = rndGen_.sample01<scalar>();

    do
    {
        rPSIm = rndGen_.sample01<scalar>();
        prob =
            pow(h1, h1)
           /(pow(h2, h2)*pow(h3, h3))
           *pow(rPSIm, h2)
           *pow(1.0 - rPSIm, h3);
    } while (prob < eps);

    return rPSIm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcCloud::dsmcCloud
(
    const Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    CloudWithModels<dsmcParcel>(mesh, cloudName, false),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(particleProperties_.get<scalar>("nEquivalentParticles")),
    axisymmetric_(particleProperties_.get<Switch>("axisymmetricSimulation")),
    radialExtent_(0.0),
    maxRWF_(1.0),
    nTerminalOutputs_
    (
        mesh.time().controlDict().get<label>("nTerminalOutputs")
    ),
    rhoNMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    momentumMeanElectron_(mesh_.nCells(), Zero),
    momentumMean_(mesh_.nCells(), Zero),
    linearKEMeanElectron_(mesh_.nCells(), 0.0),
    electronTemperature_(mesh_.nCells(), 0.0),
    cellVelocity_(mesh_.nCells(), Zero),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_(mesh_.nCells(), 0),
    constProps_(),
    rndGen_(label(clock::getTime()) + 7183*Pstream::myProcNo()),
    controllers_(t, mesh, *this),
    fields_(t, mesh, *this),
    boundaries_(t, mesh, *this),
    trackingInfo_(mesh, *this, true),
    binaryCollisionModel_
    (
        BinaryCollisionModel::New
        (
            particleProperties_,
            *this
        )
    ),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh, *this),
    boundaryMeas_(mesh, *this, true),
    cellMeas_(mesh, *this, true),
    functions_
    (
        *this,
        particleProperties_.subOrEmptyDict("cloudFunctions"),
        solution_.active()
    )
{
    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildConstProps();

    if (axisymmetric_)
    {
        radialExtent_ = particleProperties_.get<scalar>("radialExtentOfDomain");
        maxRWF_ = particleProperties_.get<scalar>("maxRadialWeightingFactor");
    }

    reactions_.initialConfiguration();

    buildCellOccupancy();

    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.sample01<scalar>();
    }

    collisionPartnerSelectionModel_ =
        collisionPartnerSelection::New(mesh, *this, particleProperties_);

    collisionPartnerSelectionModel_->initialConfiguration();

    fields_.createFields();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();
}


Foam::dsmcCloud::dsmcCloud
(
    const Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict,
    const bool clearFields
)
:
    CloudWithModels<dsmcParcel>(mesh, cloudName, false),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(particleProperties_.get<scalar>("nEquivalentParticles")),
    axisymmetric_(particleProperties_.get<Switch>("axisymmetricSimulation")),
    radialExtent_(0.0),
    maxRWF_(1.0),
    nTerminalOutputs_(mesh.time().controlDict().get<label>("nTerminalOutputs")),
    rhoNMeanElectron_(),
    rhoMMeanElectron_(),
    rhoMMean_(),
    momentumMeanElectron_(),
    momentumMean_(),
    linearKEMeanElectron_(),
    electronTemperature_(),
    cellVelocity_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_(),
    constProps_(),
    rndGen_(label(clock::getTime()) + 1526*Pstream::myProcNo()),
    controllers_(t, mesh),
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    boundaryMeas_(mesh, *this),
    cellMeas_(mesh, *this),
    functions_
    (
        *this,
        particleProperties_.subOrEmptyDict("cloudFunctions"),
        solution_.active()
    )
{
    if (!clearFields)
    {
        dsmcParcel::readFields(*this);
    }

    label initialParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(initialParcels, sumOp<label>());
    }

    if (clearFields)
    {
        Info << "clearing existing field of parcels " << endl;

        clear();

        initialParcels = 0;
    }

    if (axisymmetric_)
    {
        radialExtent_ =
            particleProperties_.get<scalar>("radialExtentOfDomain");
        maxRWF_ = particleProperties_.get<scalar>("maxRadialWeightingFactor");
    }

    buildConstProps();
    dsmcAllConfigurations conf(dsmcInitialiseDict, *this);
    conf.setInitialConfig();

    label finalParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(finalParcels, sumOp<label>());
    }

    Info<< nl << "Initial no. of parcels: " << initialParcels
        << " added parcels: " << finalParcels - initialParcels
        << ", total no. of parcels: " << finalParcels
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::dsmcCloud::getTypeIDs(const dictionary& dict) const
{
    const wordHashSet names(dict.lookup("typeIds"));

    if (names.empty())
    {
        FatalErrorInFunction
            << "Entry typeIds cannot be empty in " << dict.name()
            << exit(FatalError);
    }

    labelList typeIDs(names.size(), -1);

    label i = 0;
    forAllIters(names, iter)
    {
        const word& name = iter();

        label id = typeIdList_.find(name);

        if (id == -1)
        {
            FatalErrorInFunction
                << "Cannot find particle type: " << name << nl << "in: "
                << particleProperties_.name()
                << exit(FatalError);
        }

        typeIDs[i++] = id;
    }

    return typeIDs;
}


void Foam::dsmcCloud::evolve()
{
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);

    functions_.preEvolve();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();

    // Move the particles ballistically with their current velocities
    Cloud<dsmcParcel>::move(*this, td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    if (axisymmetric_)
    {
        axisymmetricWeighting();
        buildCellOccupancy();
    }


    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();

    // Calculate new velocities via stochastic collisions
    collisions();

    // Update cell occupancy (reactions may have changed it)
    buildCellOccupancy();

    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();

    reactions_.outputData();

    fields_.calculateFields();
    fields_.writeFields();

    controllers_.calculateProps();
    controllers_.outputResults();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    trackingInfo_.clean();
    boundaryMeas_.clean();
    cellMeas_.clean();

    functions_.postEvolve();
}


Foam::label Foam::dsmcCloud::nTerminalOutputs()
{
    return nTerminalOutputs_;
}


void Foam::dsmcCloud::info() const
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    scalar nMol = nDsmcParticles*nParticle_;

    scalar linearKineticEnergy = infoMeasurements()[1];
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar rotationalEnergy = infoMeasurements()[2];
    reduce(rotationalEnergy, sumOp<scalar>());

    scalar vibrationalEnergy = infoMeasurements()[3];
    reduce(vibrationalEnergy, sumOp<scalar>());

    scalar electronicEnergy = infoMeasurements()[4];
    reduce(electronicEnergy, sumOp<scalar>());

    Info<< "    Number of dsmc particles        = " << nDsmcParticles << endl;

    if (nDsmcParticles)
    {
        scalar totalEnergy =
            rotationalEnergy
          + linearKineticEnergy
          + vibrationalEnergy
          + electronicEnergy;

        Info<< "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average rotational energy       = "
            << rotationalEnergy/nMol << nl
            << "    Average vibrational energy      = "
            << vibrationalEnergy/nMol << nl
            << "    Average electronic energy       = "
            << electronicEnergy/nMol << nl
            << "    Total energy                    = "
            << totalEnergy
            << endl;
    }
}


Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    const scalar temperature,
    const scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.GaussNormal<vector>();
}

Foam::vector Foam::dsmcCloud::chapmanEnskogLinearVelocity
(   
    scalar breakdownParameter,
    scalar pressure,
    scalar temperature,
    scalar mass,
    vector heatFlux,
    tensor stress
)
{

    scalar mostProbableSpeed(
        maxwellianMostProbableSpeed
        (
            temperature,
            mass
        )
    );

    vector U=vector::zero;
    scalar amplitudeParameter = 1.0 + 30.0*breakdownParameter;
    scalar gamma=0.0;
    
    do {

        U=rndGen_.GaussNormal<vector>()/sqrt(2.0);

        gamma=1.0+(2.0/mostProbableSpeed*(heatFlux.x()*U.x()+heatFlux.y()*U.y()+heatFlux.z()*U.z())*(2.0/5.0*(pow(U.x(),2)+pow(U.y(),2)+pow(U.z(),2))-1.0)
             -2.0*(stress.xy()*U.x()*U.y()+stress.xz()*U.x()*U.z()+stress.yz()*U.y()*U.z())
             -stress.xx()*(pow(U.x(),2)-pow(U.z(),2))-stress.yy()*(pow(U.y(),2)-pow(U.z(),2)))/pressure;

    } while(amplitudeParameter*rndGen_.sample01<scalar>()>gamma);

    return mostProbableSpeed*U;

}

Foam::scalar Foam::dsmcCloud::equipartitionRotationalEnergy
(
    const scalar temperature,
    const scalar rotationalDof
)
{
    scalar ERot = 0.0;

    if (rotationalDof < SMALL)
    {
        return ERot;
    }
    else if (rotationalDof < 2.0 + SMALL && rotationalDof > 2.0 - SMALL)
    {
        // Special case for rDof = 2, i.e. diatomics;
        ERot =
           -log(rndGen_.sample01<scalar>())
           *physicoChemical::k.value()
           *temperature;
    }
    else
    {
        const scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        const scalar eps = rndGen_.sample01<scalar>();

        do
        {
            energyRatio = 10*rndGen_.sample01<scalar>();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < eps);

        ERot = energyRatio*physicoChemical::k.value()*temperature;
    }

    return ERot;
}


Foam::labelList Foam::dsmcCloud::equipartitionVibrationalEnergyLevel
(
    const scalar temperature,
    const scalar vibrationalDof,
    const label typeId
)
{
    labelList vibLevel(vibrationalDof, 0);

    if (vibrationalDof < SMALL)
    {
        return vibLevel;
    }
    else
    {
        forAll(vibLevel, i)
        {
            vibLevel[i] =
               -log(rndGen_.sample01<scalar>())
               *temperature
               /constProps(typeId).thetaV()[i];
        }
    }

    return vibLevel;
}


Foam::label Foam::dsmcCloud::equipartitionElectronicLevel
(
    const scalar temperature,
    const List<label>& degeneracyList,
    const List<scalar>& electronicEnergyList,
    const label typeId
)

{
    const scalar EMax = physicoChemical::k.value()*temperature;
    const label jMax = constProps(typeId).nElectronicLevels();

    label jDash = 0;

    if (jMax == 1)
    {
        return 0;
    }
    if (temperature < VSMALL)
    {
        return jDash;
    }
    else
    {
        // Calculate summation term in denominator of Eq.3.1.1
        // from Liechty thesis.
        label i = 0;
        scalar expSum = 0.0;
        do
        {
            expSum += degeneracyList[i]*exp((-electronicEnergyList[i]/EMax));
            i += 1;
        } while (i < jMax);

        // select maximum integer energy level based on boltz value.
        // Note that this depends on the temperature.

        scalar boltzMax = 0.0;

        label jSelect = 0;
        for (label ii = 0; ii < jMax; ii++)
        {
            // Eq. 3.1.1 of Liechty thesis.
            scalar boltz =
                degeneracyList[ii]
               *exp((-electronicEnergyList[ii]/EMax))
               /expSum;

            if (boltzMax < boltz)
            {
                boltzMax = boltz;
                jSelect = ii;
            }
        }

        const scalar EJ = electronicEnergyList[jSelect];
        const scalar gJ = degeneracyList[jSelect];
        const scalar expMax = gJ*exp((-EJ/EMax));

        scalar func = 0.0;
        const scalar eps = rndGen_.sample01<scalar>();

        do // acceptance-rejection based on Eq. 3.1.2 of Liechty thesis.
        {
            jDash = rndGen_.position<label>(0, jMax - 1);
            func =
                degeneracyList[jDash]
               *exp((-electronicEnergyList[jDash]/EMax))
               /expMax;
        } while (!(func > eps));
    }

    return jDash;
}


Foam::scalar Foam::dsmcCloud::postCollisionRotationalEnergy
(
    const scalar rotationalDof,
    const scalar ChiB
)
{
    scalar energyRatio = 0.0;

    if (rotationalDof == 2.0)
    {
        energyRatio = 1.0 - pow(rndGen_.sample01<scalar>(), (1.0/ChiB));
    }
    else
    {
        const scalar ChiA = 0.5*rotationalDof;
        const scalar ChiAMinusOne = ChiA - 1;
        const scalar ChiBMinusOne = ChiB - 1;

        if (ChiAMinusOne<SMALL && ChiBMinusOne<SMALL)
        {
            return rndGen_.sample01<scalar>();
        }

        scalar P;

        const scalar eps = rndGen_.sample01<scalar>();

        do
        {
            P = 0;

            energyRatio = rndGen_.sample01<scalar>();

            if (ChiAMinusOne<SMALL)
            {
                P = pow((1.0 - energyRatio), ChiBMinusOne);
            }
            else if (ChiBMinusOne<SMALL)
            {
                P = pow((1.0 - energyRatio), ChiAMinusOne);
            }
            else
            {
                P =
                    pow
                    (
                        (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                        ChiAMinusOne
                    )
                *pow
                    (
                        (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                        /ChiBMinusOne,
                        ChiBMinusOne
                    );
            }
        } while (P < eps);
    }

    return energyRatio;
}


Foam::label Foam::dsmcCloud::postCollisionVibrationalEnergyLevel
(
    const bool postReaction,
    const label vibLevel,
    const label iMax,
    const scalar thetaV,
    const scalar thetaD,
    const scalar refTempZv,
    const scalar omega,
    const scalar Zref,
    const scalar Ec
)
{
    label iDash = vibLevel;

    if (postReaction)
    {
        // post - collision quantum number
        scalar func = 0.0;
        scalar EVib = 0.0;

        do // acceptance-rejection
        {
            iDash = rndGen_.position<label>(0, iMax);
            EVib = iDash*physicoChemical::k.value()*thetaV;

            // equation 5.61, Bird
            func = pow((1.0 - (EVib / Ec)), (1.5 - omega));

        } while (!(func > rndGen_.sample01<scalar>()));
    }
    else
    {
        // "quantised collision temperature" (equation 3, Bird 2010),
        // denominator from Bird 5.42

        const scalar TColl = (iMax*thetaV)/(3.5 - omega);
        const scalar pow1 = pow((thetaD/TColl), 0.33333) - 1.0;
        const scalar pow2 = pow((thetaD/refTempZv), 0.33333) - 1.0;

        // vibrational collision number (equation 2, Bird 2010)
        const scalar ZvP1 = pow((thetaD/TColl), omega);
        const scalar ZvP2 =
            pow(Zref*pow(thetaD/refTempZv, -omega), pow1/pow2);
        const scalar Zv = ZvP1*ZvP2;

        // In order to obtain the relaxation rate corresponding to Zv with the
        // collision energy - based procedure, the inelastic fraction should be
        // set to about 1/(5Zv) Bird 2008 RGD "A Comparison of Collison
        // Energy - Based and Temperature-Based..."

        const scalar inverseVibrationalCollisionNumber = 1.0/(5.0*Zv);

        if (inverseVibrationalCollisionNumber > rndGen_.sample01<scalar>())
        {
            // post - collision quantum number
            scalar func = 0.0;
            scalar EVib = 0.0;

            do // acceptance-rejection
            {
                iDash = rndGen_.position<label>(0, iMax);
                EVib = iDash*physicoChemical::k.value()*thetaV;

                // equation 5.61, Bird
                func = pow((1.0 - (EVib / Ec)), (1.5 - omega));

            } while (!(func > rndGen_.sample01<scalar>()));
        }
    }

    return iDash;
}


Foam::label Foam::dsmcCloud::postCollisionElectronicEnergyLevel
(
    const scalar Ec,
    const label jMax,
    const scalar omega,
    const List<scalar>& EElist,
    const List<label>& gList
)
{
    label nPossStates = 0;
    label ELevel = -1;

    if (jMax == 1)
    {
        nPossStates = gList[0];
    }
    else
    {
        forAll(EElist, n)
        {
            if (Ec > EElist[n])
            {
                nPossStates += gList[n];
            }
        }
    }

    label II = 0;

    do
    {
        const label nState = ceil(rndGen_.sample01<scalar>()*(nPossStates));
        label nAvailableStates = 0;
        label nLevel = -1;

        forAll(EElist, n)
        {
            nAvailableStates += gList[n];

            if (nState <= nAvailableStates && nLevel < 0)
            {
                nLevel = n;
            }
        }

        if (Ec > EElist[nLevel])
        {
            scalar prob = pow(1.0 - (EElist[nLevel]/Ec), 1.5 - omega);

            if (prob > rndGen_.sample01<scalar>())
            {
                II = 1;
                ELevel = nLevel;
            }
        }

    } while (II == 0);

    return ELevel;
}


void Foam::dsmcCloud::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        const vector position(p.position());

        pObj<< "v " << position.x()
            << " "  << position.y()
            << " "  << position.z()
            << nl;
    }

    pObj.flush();
}


void Foam::dsmcCloud::axisymmetricWeighting()
{
    const auto& cellOccupancy = this->cellOccupancy();

    forAll(cellOccupancy, celli)
    {
        for (dsmcParcel* p : cellOccupancy[celli])
        {
            const point& cC = mesh_.cellCentres()[celli];
            const scalar radius = cC.y();
            const scalar oldRadialWeight = p->RWF();
            const scalar newRadialWeight = 1.0 + maxRWF_*(radius/radialExtent_);

            p->RWF() = newRadialWeight;

            if (oldRadialWeight > newRadialWeight)
            {
                // particle might be cloned

                scalar prob = (oldRadialWeight/newRadialWeight) - 1.0;

                while (prob > 1.0)
                {
                    // add a particle and reduce prob by 1.0

                    const vector position(p->position());

                    vector U = p->U();

                    U.z() *= -1.0;

                    addNewParcel
                    (
                        position,
                        U,
                        p->RWF(),
                        p->ERot(),
                        p->ELevel(),
                        celli,
                        p->typeId(),
                        p->newParcel(),
                        p->vibLevel()
                    );

                    prob -= 1.0;
                }

                if (prob > rndGen_.sample01<scalar>())
                {
                    const vector position(p->position());

                    vector U = p->U();

                    U.z() *= -1.0;

                    addNewParcel
                    (
                        position,
                        U,
                        p->RWF(),
                        p->ERot(),
                        p->ELevel(),
                        celli,
                        p->typeId(),
                        p->newParcel(),
                        p->vibLevel()
                    );
                }
            }

            if (newRadialWeight > oldRadialWeight)
            {
                // particle might be deleted
                if
                (
                    oldRadialWeight/newRadialWeight
                  < rndGen_.sample01<scalar>()
                )
                {
                    deleteParticle(*p);
                }
            }
        }
    }
}


// ************************************************************************* //
