/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "xParticleCloud.H"
// #include "zeroGradientFvPatchFields.H"
#include "constants.H"
// #include "GeometricField.H"
#include "fvMesh.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::xParticleCloud::buildConstProps()
{
    Info<< nl <<
    "Constructing constant properties for eatraterrestrial particles." << endl;

    constProps_.setSize(typeIdList_.size());

    dictionary particleProperties
    (
        particleProperties_.subDict("particleProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& particleDict(particleProperties.subDict(id));

        constProps_[i] = xParticle::constantProperties(particleDict);
    }
}

void Foam::xParticleCloud::moveCollide
(
    xParticleCloud& cloud,
    xParticle::trackingData& td,
    const scalar deltaT
)
{
    xParticle::trackingData td0(*this, 0);
    Cloud<xParticle>::move(*this, td0, mesh_.time().deltaTValue());

    xParticle::trackingData td1(*this, 1);
    Cloud<xParticle>::move(*this, td0, mesh_.time().deltaTValue());

    this->updateCellOccupancy();

    //- update particle force
    calculateForce();

    //- collision() returns collisionModel
    this->collision().collide();

    // xParticle::trackingData td0(*this, 0);
    Cloud<xParticle>::move(*this, td0, mesh_.time().deltaTValue());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//-Constructer of running simulations
Foam::xParticleCloud::xParticleCloud
(
    const Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    dsmcCloud* dsmcCloudPtr,
    bool readFields
)
:
    CloudWithModels<xParticle>(mesh, cloudName, false),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nRealParticles_(particleProperties_.get<scalar>("nRealParticles")),
    maxDistance_(particleProperties_.get<scalar>("maxPotentialContactDistance")),
    constProps_(),
    il_(mesh,maxDistance_,false),
    collisionModel_(nullptr),
    dsmcCloudPtr_(dsmcCloudPtr),
    rndGen_(label(clock::getTime()) + 7183*Pstream::myProcNo())
{   
    if (readFields)
    {
        xParticle::readFields(*this);
    }
    
    buildConstProps();
    
    buildCellOccupancy();

    //- setup collision model
    collisionModel_.reset
    (
        //- can be found in agrangian/intermediate/submodels/Kinematic/CollisionModel/CollisionModel/CollisionModelNew.C

        CollisionModel<xParticleCloud>::New
        (
            particleProperties_.subDict("interparticleCollisions"),
            *this
        ).ptr()
    );
}

//- Constructer of Initialization
Foam::xParticleCloud::xParticleCloud
(
    const Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& initialiseDict,
    dsmcCloud* dsmcCloudPtr,
    const bool& clearFields
)
    :
    CloudWithModels<xParticle>(mesh, cloudName, false),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nRealParticles_(particleProperties_.get<scalar>("nRealParticles")),
    maxDistance_(particleProperties_.get<scalar>("maxPotentialContactDistance")),
    constProps_(),
    il_(mesh,0,false),
    collisionModel_(nullptr),
    dsmcCloudPtr_(dsmcCloudPtr),
    rndGen_(label(clock::getTime()) + 1526*Pstream::myProcNo())
{
    if(!clearFields)
    {
        xParticle::readFields(*this);
    }

//     label initialSolidParticles = this->size();

//     if (Pstream::parRun())
//     {
//         reduce(initialSolidParticles, sumOp<label>());
//     }
    
    if(clearFields)
    {
        Info << "Clearing existing field of solid particles " << endl;

        clearParticles();

//         initialSolidParticles = 0;
    }
    
    buildConstProps();

//     solidAllConfigurations conf(solidInitialiseDict, *this);
//     conf.setInitialConfig();
    
    
//     if (Pstream::parRun())
//     {
//         reduce(finalSolidParticles, sumOp<label>());
//     }

//     Info << nl << "Initial no. of solid particles: " << initialSolidParticles
//          << " added solid particles: " << finalSolidParticles - initialSolidParticles
//          << ", total no. of solid particles: " << finalSolidParticles
//          << endl;

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::labelList Foam::xParticleCloud::getTypeIDs(const dictionary& dict) const
{
    const wordHashSet names(dict.lookup("typeIds"));

    if (names.empty())
    {
        FatalErrorInFunction
            << "Entry solid phase typeIds cannot be empty in " << dict.name()
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

/*---------- Begin of Evolve ----------*/
// -- Evolve function used in rarefiedMultiphaseFoam.C
void Foam::xParticleCloud::evolve()
{

    motion();
    info();
    
}
/*---------- End of Evolve ----------*/

void Foam::xParticleCloud::calculateForce()
{

}

void Foam::xParticleCloud::motion
(
    xParticleCloud& cloud,
    xParticle::trackingData& td
)
{
    label nSubCycles = collision().nSubCycles();

    if (nSubCycles > 1)
    {
        Info<< "    " << nSubCycles << " move-collide subCycles" << endl;

        subCycleTime moveCollideSubCycle
        (
            const_cast<Time&>(this->db().time()),
            nSubCycles
        );

        while(!(++moveCollideSubCycle).end())
        {
            moveCollide(cloud, td, this->db().time().deltaTValue());
        }

        moveCollideSubCycle.endSubCycle();
    }
    else
    {
        moveCollide(cloud, td, this->db().time().deltaTValue());
    }
}

void Foam::xParticleCloud::axisymmetricWeighting()
{
    const auto& cellOccupancy = this->cellOccupancy();

    forAll(cellOccupancy, c)
    {
        for (xParticle* p : cellOccupancy[c])
        {
            const point& cC = mesh_.cellCentres()[c];
            const scalar radius = cC.y();
                        
            const scalar oldRadialWeight = p->RWF();
            const scalar newRadialWeight = 1.0 + dsmcCloudPtr()->maxRWF()*(radius/dsmcCloudPtr()->radialExtent());

            p->RWF() = newRadialWeight;
                
            if(oldRadialWeight > newRadialWeight) 
            {
                //particle might be cloned
                    
                scalar prob = (oldRadialWeight/newRadialWeight) - 1.0;
                    
                while(prob > 1.0)
                {
                    //add a particle and reduce prob by 1.0
                    
                    const vector position(p->position());
                    
                    vector Us = p->U();
                    
                    Us.z() *= -1.0;

                    addNewParticle
                    (
                        position,
                        c,
                        p->typeID(),
                        p->newParticle(),
                        p->T(),
                        p->RWF(),
                        Us,
                        // p->omega(),
                        p->f(),
                        p->angularMomentum(),
                        p->torque()
                        
                    );
                        
                    prob -= 1.0;
                }
                    
                if(prob > rndGen_.sample01<scalar>())
                {
                    const vector position = p->position();
                    
                    vector Us = p->U();
                    
                    Us.z() *= -1.0;

                    addNewParticle
                    (
                        position,
                        c,
                        p->typeID(),
                        p->newParticle(),
                        p->T(),
                        p->RWF(),
                        Us,
                        // p->omega(),
                        p->f(),
                        p->angularMomentum(),
                        p->torque()
                    );
                }
            }
            
            if(newRadialWeight > oldRadialWeight)
            {           
                //particle might be deleted
                if((oldRadialWeight/newRadialWeight) < rndGen_.sample01<scalar>())
                {
                    deleteParticle(*p);
                } 
            } 
        }
    }
}

void Foam::xParticleCloud::addNewParticle
(
    const vector& position,
    const label cellI,
    const label typeID,
    const label newParticle,
    const scalar T,
    const scalar RWF,
    const vector& U,
    // const vector& omega,
    const vector& f,
    const vector& angularMomentum,
    const vector& torque
)
{
    xParticle* pPtr = new xParticle
    (
        mesh_,
        position,
        cellI,
        typeID,
        newParticle,
        T,
        RWF,
        U,
        // omega,
        f,
        angularMomentum,
        torque
    );

    addParticle(pPtr);
}

void Foam::xParticleCloud::info() const
{
    label nSimulators = this->size();
    reduce(nSimulators, sumOp<label>());

    
    Info << "    Number of xt particles       = " << nSimulators<< endl;
    
    Info << endl;
    
}

// ************************************************************************* //
