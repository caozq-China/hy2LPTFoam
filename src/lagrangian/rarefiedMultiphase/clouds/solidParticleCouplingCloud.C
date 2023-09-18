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

#include "solidParticleCouplingCloud.H"
#include "zeroGradientFvPatchFields.H"
#include "constants.H"
// #include "GeometricField.H"

using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<solidParticleCoupling>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidParticleCouplingCloud::enableMoleculeVelocityUpdate()
{
    return interphaseCouplingModel_->enableMoleculeVelocityUpdate();
}

// bool Foam::solidParticleCouplingCloud::enableStochasticCollision()
// {
//     return solidCollisionDetectionModel_->active();
// }

void Foam::solidParticleCouplingCloud::interphaseCouplingCal()
{
    label nInterphaseCollisions = 0;
    forAll(cellOccupancy(), cellI)
    {
        const DynamicList<solidParticleCoupling*>& cellsolidParcels =  cellOccupancy()[cellI];
        const DynamicList<dsmcParcel*>& cellParcels = dsmcCloudReference_->cellOccupancy()[cellI];
        //- if there are gas and solid particles in a cell at the same time,
        //- then the interphase effect must be considered.
        if(cellsolidParcels.size() > 0 && cellParcels.size() > 0)
        {
            
            scalar gasRWF = 0.0;
            
            //- calculating RWF for gas phase
            if(dsmcCloudReference_->axisymmetric())
            {
                scalar nMols = 0.0;
                
                forAll(cellParcels, i)
                {
                    const dsmcParcel& pMol = *cellParcels[i];
                    
                    scalar radiusGas = sqrt(sqr(pMol.position().y()) 
                                        + sqr(pMol.position().z()));

                    gasRWF += 1.0 + dsmcCloudReference_->maxRWF()*(radiusGas/dsmcCloudReference_->radialExtent());
                    
                    nMols += 1.0;
                }
                
                gasRWF /= nMols;
            }
            else
            {
                gasRWF=1.0;
            }
            
            //- reset momentum transfer
//             forAll(cellsolidParcels,solidParticleID)
//             {
//                 solidParticleCoupling& pSolid = *cellsolidParcels[solidParticleID];
//                 pSolid.F() = Zero;
//             }
        
            //- loop over all solid particles in the cell
            forAll(cellsolidParcels,solidParticleID)
            {
                //- solid particle properties
                solidParticleCoupling& pSolid = *cellsolidParcels[solidParticleID];
                const scalar& cellVolume = mesh_.cellVolumes()[cellI];
                
                //- time step
                const scalar deltaT = mesh_.time().deltaTValue();
                
                //- reset momentum transfer
                pSolid.F() = Zero;
                
                // ------ 1. CALCULATE THE MOMENTUM AND ENERGY TRANSFER ------ //
                vector FdeltaC(Zero);
                scalar QdeltaC = 0.0;
                    
                //- loop over all dsmcparcels in the cell to
                //- calculate enrgy and force addition
                forAll(cellParcels, i)
                {
                    dsmcParcel& p = *cellParcels[i];
                    const label typeIdP = p.typeId();
                    const scalar rotationalDof = dsmcCloudReference_->constProps(typeIdP).rotationalDoF();
                    
                    if(rotationalDof == 0)
                    {
                        //- if it is the case of gas atoms
                        FdeltaC += interphaseCoupling().monatomicForceTransferToSolidParticle
                        (
                            pSolid,
                            p,
                            cellVolume,
                            gasRWF
                        );
                        
                        QdeltaC += interphaseCoupling().monatomicEnergyTransferToSolidParticle
                        (
                            pSolid,
                            p,
                            cellVolume,
                            gasRWF
                        );
                    }
                    else
                    {
                        //- if it is the case of gas molecules
                        FdeltaC += interphaseCoupling().polyatomicForceTransferToSolidParticle
                        (
                            pSolid,
                            p,
                            cellVolume,
                            gasRWF
                        );
                        
                        QdeltaC += interphaseCoupling().polyatomicEnergyTransferToSolidParticle
                        (
                            pSolid,
                            p,
                            cellVolume,
                            gasRWF
                        );
                    }
                }
                
                //Info<<"FdeltaC = "<<FdeltaC.x() << endl;
            
                /* UPDATING SOLID PARTICLE PROPERTIES */
                pSolid.F() = FdeltaC;
                
                solidPhaseChange().temperatureCorrection(pSolid, QdeltaC, deltaT);
                
                // End of Momentum and Energy Transfer
                
                
                // ------ 2. CALCULATE THE GAS ATOM/MOLECULE REFLECTION ------//
                if(enableMoleculeVelocityUpdate())
                {
                    
                    scalar sigmaTcRMax = sigmaTcRMax_[cellI];
                    //- reset the number of reflection candidae
                    nDsmcReflectionCandidate_ = 0.0;
                    
                    //scalar dsolid = constSolidProps(typeIdSolid).dSolid();
                    const scalar dsolid = pSolid.D();
                    
                    if(dsmcCloudReference_->axisymmetric())
                    {
                        scalar solidRWF = 0.0;
                        scalar nParticles = 0.0;
                        
                        forAll(cellsolidParcels, i)
                        {
                            const solidParticleCoupling& pSolid = *cellsolidParcels[i];
                            
                            scalar radiusSolid = sqrt(sqr(pSolid.position().y()) 
                                                + sqr(pSolid.position().z()));

                            solidRWF += 1.0 + dsmcCloudReference_->maxRWF()*(radiusSolid/dsmcCloudReference_->radialExtent());
                            
                            nParticles += 1.0;
                        }
                        
                        solidRWF /= nParticles;
                        
                        nDsmcReflectionCandidate_ = nSolidParticles_*solidRWF*cellParcels.size()*
                                sigmaTcRMax*deltaT/cellVolume;
                        
                    }
                    else
                    {
                        nDsmcReflectionCandidate_ = nSolidParticles_*cellParcels.size()*
                                sigmaTcRMax*deltaT/cellVolume;
                    }
                    
                    //- label always round down nDsmcReflectionCandidate_ so a random number between 0 and 1 is added
                    //- to give an even possibility to round up
                    label nDsmcCandidate(nDsmcReflectionCandidate_+rndGenS_.sample01<scalar>());
                    
                    if(nDsmcCandidate > 0)
                    {
                        labelList dsmcParcelLabel(cellParcels.size());
                        for(label i=0; i < cellParcels.size(); ++i)
                        {
                            dsmcParcelLabel[i] = i;
                        }
                        /*--------------------------------------------------------------------*\
                        cellParcels: [a,b,c,d,e,f,g...]
                        dsmcParcelLabel: [0,1,2,3,4,5,6,7,8...]
                        when one is picked, such as 4 in dsmcParcelLabel 
                        and the corresponding parcel is number dsmcParcelLabel[4] which is e
                        \*--------------------------------------------------------------------*/
                        
                        label selectedDsmcCandidate = -1;
                        
                        for(label i=0; i < nDsmcCandidate; ++i)
                        {
                            //- pick a dsmcparcel in the cellParcels
                            selectedDsmcCandidate = rndGenS_.position<label>(0, (cellParcels.size()-1));
                            
                            dsmcParcel& potentialDsmcCandidate = *cellParcels[dsmcParcelLabel[selectedDsmcCandidate]];
                            
                            
                            scalar sigmaTcRp = pi*sqr(dsolid/2.0)*mag(potentialDsmcCandidate.U()-pSolid.U());
                            
                            //Info<<"cRp = "<< cRp<<endl;
                            
                            if(sigmaTcRp > sigmaTcRMax)
                            {
                                sigmaTcRMax_[cellI] = sigmaTcRp;
                            }
                            
                            if((sigmaTcRp/sigmaTcRMax) > rndGenS_.sample01<scalar>())
                            {
                                //Info<<"Doing interphase collide"<<endl;
                                interphaseCoupling().moleculePostCollisionVelocityUpdate
                                (
                                    pSolid,
                                    potentialDsmcCandidate
                                );
                                
                                solidCellMeas_.nColls()[cellI]++;
                                solidCellMeas_.nCollsTotal()[Pstream::myProcNo()]++;
                                
                                ++nInterphaseCollisions;
                                /*
                                //- To get 1 million samples (for validation only)
                                if(nInterphaseCollisions==1000000)
                                {
                                    //- force to stop the simulation
                                    FatalErrorInFunction
                                    << "    Sampling Done! "
                                    << nl
                                    << exit(FatalError); 
                                }
                                */
                                
                            }
                        }
                    }
                }// End of Gas Reflection
            }
        }
    }
    
    reduce(nInterphaseCollisions, sumOp<label>());
    
    Info << "    Interphase collisions           = " << nInterphaseCollisions << endl;
    sigmaTcRMax_.correctBoundaryConditions();
}


void Foam::solidParticleCouplingCloud::buildSolidConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constSolidProps_.setSize(typeIdSolidList_.size());

    dictionary solidProperties
    (
        particleProperties_.subDict("solidProperties")
    );

    forAll(typeIdSolidList_, i)
    {
        const word& id(typeIdSolidList_[i]);

        Info<< "    " << id << endl;

        const dictionary& solidDict(solidProperties.subDict(id));

        constSolidProps_[i] = solidParticleCoupling::constantProperties(solidDict);
    }
}

//- beginning of solid particle velocity update -//
void Foam::solidParticleCouplingCloud::updateVelocityHalfTimeStep()
{
    
    forAllIter(solidParticleCouplingCloud, *this, iter)
    {
        solidParticleCoupling* p = &iter();
        
        p->U() += 0.5*(p->F()/constSolidProps(p->typeID()).massSphere())*mesh_.time().deltaTValue();
    }
}

void Foam::solidParticleCouplingCloud::updateVelocityHalfTimeStepBack()
{
    
    forAllIter(solidParticleCouplingCloud, *this, iter)
    {
        solidParticleCoupling* p = &iter();
        
        p->U() -= 0.5*(p->F()/constSolidProps(p->typeID()).massSphere())*mesh_.time().deltaTValue();
    }
}

void Foam::solidParticleCouplingCloud::updateVelocityFullTimeStep()
{
    
    forAllIter(solidParticleCouplingCloud, *this, iter)
    {
        solidParticleCoupling* p = &iter();
        
        p->U() += (p->F()/constSolidProps(p->typeID()).massSphere())*mesh_.time().deltaTValue();
    }
}

void Foam::solidParticleCouplingCloud::updateMppicCellList()
{
    forAll(cellOccupancy(),cellI)
    {
        
        if(volumeFraction_[cellI] > criticalVolFrac_)
        {
            cellSolidPackingOccupIds_.append(cellI);
        }
        else
        {
            cellSolidCollisionOccupIds_.append(cellI);
        }
    }
}


void Foam::solidParticleCouplingCloud::clearMppicCellList()
{
    cellSolidPackingOccupIds_.clear();
//     cellSolidPackingOccupIds_.shrink();
    cellSolidCollisionOccupIds_.clear();
//     cellSolidCollisionOccupIds_.shrink();
}

void Foam::solidParticleCouplingCloud::particleSizeCorrection
(
    solidParticleCoupling& pSolid
)
{
    
    if(pSolid.phaseState() == 2)
    {
        //- the core is solid
        scalar rhoLiquid = particleLiquidDensityCorrection
                            (
                                materialList_[pSolid.typeID()],
                                pSolid.T()
                            );
        pSolid.D() = pow(
                            (6.0*constSolidProps(pSolid.typeID()).massSphere()/pi)
                            /((1.0-pow(pSolid.CzRatio(),3.0))*constSolidProps(pSolid.typeID()).rho()
                            +pow(pSolid.CzRatio(),3.0)*rhoLiquid)
                            ,1.0/3.0
                         );
    }
    else if(pSolid.phaseState() == 3 && pSolid.T() >= constSolidProps(pSolid.typeID()).Tm())
    {
        //- pure liquid droplet
        scalar rhoLiquid = particleLiquidDensityCorrection
                            (
                                materialList_[pSolid.typeID()],
                                pSolid.T()
                            );
        
        pSolid.D() = pow((6.0*constSolidProps(pSolid.typeID()).massSphere()/(pi*rhoLiquid)),1.0/3.0);
        
    }
    else
    {
        pSolid.D() = constSolidProps(pSolid.typeID()).d();
    }
     
}

void Foam::solidParticleCouplingCloud::updateParticleVolumeFraction()
{
    //- reset
    forAll(cellOccupancy(),cellI)
    {
        volumeFraction_[cellI] = 0.0;
    }
    //- add all particle volume
    forAllIter(solidParticleCouplingCloud, *this, iter)
    {
        solidParticleCoupling* p = &iter();
        
        scalar RWF = dsmcCloudReference()->axiRWF(mesh().cellCentres()[p->cell()]);
        
        volumeFraction_[p->cell()] += ((pi*pow(p->D(),3.0)/6.0)*RWF*nSolidParticles());
    }
    
    forAll(cellOccupancy(),cellI)
    {
        const scalar cellVolume = mesh_.cellVolumes()[cellI];
        
        volumeFraction_[cellI] /= cellVolume;
    }
}

void Foam::solidParticleCouplingCloud::particleParticleCollisions()
{
    solidCollisionDetectionModel_->collide();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//-Constructer of running simulations
Foam::solidParticleCouplingCloud::solidParticleCouplingCloud
(
    const Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    dsmcCloud* dsmcCloudReference,
    bool readFields
)
:
    CloudWithModels<solidParticleCoupling>(mesh, cloudName, false),
    mppicPropertiesDict_
    (
        IOobject
        (
            "mppicPropertiesDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdSolidList_(particleProperties_.lookup("typeIdList")),
    materialList_(typeIdSolidList_.size(),"null"),
    nSolidParticles_(particleProperties_.get<scalar>("nEquivalentSolidParticles")),
    nDsmcReflectionCandidate_(0.0),
    solidWeightFactor_(particleProperties_.get<scalar>("solidWeightingAmplificationFactor")),
    criticalVolFrac_(0.0),
    enableMPPICMethod_(mppicPropertiesDict_.get<Switch>("enableMPPICMethod")),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "sigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    cRsigmaPPMax_
    (
        IOobject
        (
            this->name() + "cRsigmaPPMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    volumeFraction_
    (
        IOobject
        (
            this->name() + "volumeFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    constSolidProps_(),
    rndGenS_(label(clock::getTime()) + 7183*Pstream::myProcNo()),
    controllers_(t, mesh, *this),
    fields_(t, mesh, *this),
    boundaries_(t, mesh, *this),
    interphaseCouplingModel_(),
    solidCollisionDetectionModel_(),
    solidBinaryCollisionModel_(),
    solidPhaseChange_(),
    cellSolidPackingOccupIds_(),
    cellSolidCollisionOccupIds_(),
    mppicDampingModel_(),
    mppicPackingModel_(),
    mppicIsotropyModel_(),
    solidBoundaryMeas_(mesh, *this, true),
    solidCellMeas_(mesh, *this, true),
    dsmcCloudReference_(dsmcCloudReference)
{   
    if (readFields)
    {
        solidParticleCoupling::readFields(*this);
    }
    
    buildSolidConstProps();
    
    solidPhaseChange_ = 
            solidPhaseChangeModel::New
            (
                *this,
                particleProperties_
            );
    
    if(solidPhaseChange_->active())
    {
        
        const List<word> MaterialLists(particleProperties_.lookup("materialList"));
        
        forAll(MaterialLists,MaterialID)
        {
            materialList_[MaterialID] = MaterialLists[MaterialID];
        }
    
    }
    
    interphaseCouplingModel_ = 
        InterphaseCoupling::New
        (
            *this,
            particleProperties_
        );
    
    buildCellOccupancy();
    
    solidCollisionDetectionModel_ = 
            solidCollisionDetection::New(mesh, *this, particleProperties_);

    solidCollisionDetectionModel_->initialConfiguration();
    
    if(solidCollisionDetectionModel_->active() == true && enableMPPICMethod_ == true)
    {
        FatalErrorInFunction
        << "    Eorror in contstructer for intialization. "
        << "    Stochastic particle-particle collision method and MPPIC method are "
        << "    enabled at the same time ! Please select one!"
        << nl
        << exit(FatalError);
    }
    else if(solidCollisionDetectionModel_->active()==true && enableMPPICMethod_ == false)
    {
        Info<<nl<<"The stochastic collision method is enabled !    "<<endl;
        solidBinaryCollisionModel_ =  
            solidBinaryCollisionModel::New
            (
                particleProperties_,
                *this
            );
    }
    else if(solidCollisionDetectionModel_->active() == false && enableMPPICMethod_ == true)
    {
        Info<<nl<<"The MPPIC method is enabled !    "<<endl;
        
        mppicDampingModel_ = 
            mppicDampingModel::New
            (
                *this,
                mppicPropertiesDict_
            );
        
        
        mppicPackingModel_ =
            mppicPackingModel::New
            (
                *this,
                mppicPropertiesDict_
            );
        
        mppicIsotropyModel_ =
            mppicIsotropyModel::New
            (
                *this,
                mppicPropertiesDict_
            );
      
        criticalVolFrac_ = particleProperties_.get<scalar>("criticalParticleVolumeFraction");
        
        word dampingModelTypeName = 
        mppicPropertiesDict_.getOrDefault<word>("mppicDampingModel", "NoMppicDamping");
        word packingModelTypeName =
        mppicPropertiesDict_.getOrDefault<word>("mppicPackingModel", "NoPacking");
        word isotropyModelTypeName =
        mppicPropertiesDict_.getOrDefault<word>("mppicIsotropyModel", "NoIsotropyModel");
//         Info<<"packingModelTypeName = "<<packingModelTypeName<<endl;
        if
        (
            dampingModelTypeName == "NoMppicDamping" && 
            packingModelTypeName == "NoPacking" &&
            isotropyModelTypeName == "NoIsotropyModel"
        )
        {
            FatalErrorInFunction
            << "    The MPPIC method is activated, but no MPPIC model is selected. "
            << "    Please check if there is any typo in system/mppicPropertiesDict "
            << nl
            << exit(FatalError);
        }
    }
    
    
    fields_.createFields();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();//- do nothing in solidGravitationalAccelerationController
}

//- Constructer of Initialization
Foam::solidParticleCouplingCloud::solidParticleCouplingCloud
(
    const Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& solidInitialiseDict,
    dsmcCloud* dsmcCloudReference, 
    const bool& clearFields
)
    :
    CloudWithModels<solidParticleCoupling>(mesh, cloudName, false),
    mppicPropertiesDict_
    (
        IOobject
        (
            "mppicPropertiesDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdSolidList_(particleProperties_.lookup("typeIdList")),
    materialList_(typeIdSolidList_.size(),"null"),
    nSolidParticles_(particleProperties_.get<scalar>("nEquivalentSolidParticles")),
    nDsmcReflectionCandidate_(0.0),
    solidWeightFactor_(particleProperties_.get<scalar>("solidWeightingAmplificationFactor")),
    criticalVolFrac_(0.0),
    enableMPPICMethod_(mppicPropertiesDict_.get<Switch>("enableMPPICMethod")),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "sigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    cRsigmaPPMax_
    (
        IOobject
        (
            this->name() + "cRsigmaPPMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    volumeFraction_
    (
        IOobject
        (
            this->name() + "volumeFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    constSolidProps_(),
    rndGenS_(label(clock::getTime()) + 1526*Pstream::myProcNo()),
    controllers_(t, mesh),
    fields_(t, mesh),
    boundaries_(t, mesh),
    interphaseCouplingModel_(),
    solidCollisionDetectionModel_(),
    solidBinaryCollisionModel_(),
    solidPhaseChange_(),
    cellSolidPackingOccupIds_(),
    cellSolidCollisionOccupIds_(),
    mppicDampingModel_(),
    mppicPackingModel_(),
    mppicIsotropyModel_(),
    solidBoundaryMeas_(mesh, *this),
    solidCellMeas_(mesh, *this),
    dsmcCloudReference_(dsmcCloudReference)
{
    if(!clearFields)
    {
        solidParticleCoupling::readFields(*this);
    }

    label initialSolidParticles = this->size();

    if (Pstream::parRun())
    {
        reduce(initialSolidParticles, sumOp<label>());
    }
    
    solidPhaseChange_ = 
            solidPhaseChangeModel::New
            (
                *this,
                particleProperties_
            );
    
    solidCollisionDetectionModel_ = 
            solidCollisionDetection::New(mesh, *this, particleProperties_);

    solidCollisionDetectionModel_->initialConfiguration();
    
    //- check if two methods of calculating particle-particle interactions are activated at the same time
    if(solidCollisionDetectionModel_->active() == true && enableMPPICMethod_ == true)
    {
        FatalErrorInFunction
        << "    Eorror in contstructer for intialization. "
        << "    Stochastic particle-particle collision method and MPPIC method are "
        << "    enabled at the same time ! Please select one!"
        << nl
        << exit(FatalError);
    }
    
    if(solidCollisionDetectionModel_->active())
    {
        Info<<"Attention:  Phase Change Model is enabled !    "<<endl;
        
        const List<word> MaterialLists(particleProperties_.lookup("materialList"));
        
        Info<<"MaterialLists = "<<MaterialLists<<endl;
        
        if(materialList_.size() != MaterialLists.size())
        {
            FatalErrorInFunction
            << "    Eorror in contstructer for intialization. "
            << "    The length of list of material does not equal to that of the list of typeId in constant/spcProperties !"
            << nl
            << exit(FatalError);
        }
        
        
        forAll(MaterialLists,MaterialID)
        {
            materialList_[MaterialID] = MaterialLists[MaterialID];
        }
        
        Info<<"MaterialList_ = "<<materialList_<<endl;
        //materialList_(particleProperties_.lookup("materialList"));
        
    }
    else
    {
        Info<<"Attention:  Phase Change Model is unabled !    "<<endl;
    }
    
    if(clearFields)
    {
        Info << "Clearing existing field of solid particles " << endl;

        clearSolidParticles();

        initialSolidParticles = 0;
    }
    
    buildSolidConstProps();
    solidAllConfigurations conf(solidInitialiseDict, *this);
    conf.setInitialConfig();
    
    label finalSolidParticles = this->size();
    
    
    if (Pstream::parRun())
    {
        reduce(finalSolidParticles, sumOp<label>());
    }

    Info << nl << "Initial no. of solid particles: " << initialSolidParticles 
         << " added solid particles: " << finalSolidParticles - initialSolidParticles
         << ", total no. of solid particles: " << finalSolidParticles 
         << endl;

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::labelList Foam::solidParticleCouplingCloud::getTypeIDs(const dictionary& dict) const
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

        label id = typeIdSolidList_.find(name);

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


Foam::scalar Foam::solidParticleCouplingCloud::conductiveEnergyTransferWithWallPatch
(
    const solidParticleCoupling& p,
    const scalar Twall,
    const scalar poissonRatioWall,
    const scalar elasticModuliWall,
    const scalar materialDensityWall,
    const scalar CpWall,
    const scalar kWall
)
{
    //- This model is based on "Heat Transfer Between Colliding Surfaces and Particles" from "Journal of Heat Transfer"
    //- by Like Li, RenWei Mei, James F.Klausner and David W.Hahn in 2012
    
//     if(enableParticlePhaseChangeModel())
//     {
//         FatalErrorInFunction
//         << "    Warning! Particle phase change model is enabled ! "
//         << "    This function is only suitable for particles in solid phase, not liquid phase. "
//         << "    Please check!"
//         << nl
//         << exit(FatalError); 
//     }
//     
//     
//     //-calculating heat transfer
//     const scalar k = constSolidProps(p.typeID()).k();
//     
//     const scalar nu = constSolidProps(p.typeID()).poissonRatio();
//     
//     const scalar E = constSolidProps(p.typeID()).elasticModuli();
//     
//     const scalar CpSolid = constSolidProps(p.typeID()).Cp();
//     
//     const scalar rhoSolid = constSolidProps(p.typeID()).rho();
//     
//     
//     const scalar massSolid = constSolidProps(p.typeID()).massSphere();
//     
//     const scalar radius = constSolidProps(p.typeID()).d()/2.0;
//     
//     const vector nw = p.normal();// -surface normal vector
//     
//     //- component in the normal direction
//     //- equal to the relative velocity between the solid particle and the wall
//     const scalar U_dot_nw = p.U() & nw;
//     
//     const scalar E12 = (4.0/3.0)/((1.0-nu*nu)/E+(1.0-sqr(poissonRatioWall))/elasticModuliWall);
//     
//     const scalar tc = 2.9432752*pow((5.0*massSolid/(4.0*E12)),0.4)*pow((radius*U_dot_nw),-0.2);
//     
//     const scalar Rc = pow((5.0*massSolid*radius*radius/(4.0*E12)),0.2)*pow(U_dot_nw,0.4);
//     
//     scalar thermalDiffusivity = k/(CpSolid*rhoSolid);
//     
//     scalar thermalDiffusivityWall = kWall/(CpWall*materialDensityWall);
//     
//     const scalar Ac = pi*Rc*Rc;
//     
//     scalar alphaPrime = thermalDiffusivityWall/thermalDiffusivity;
//     
//     scalar kPrime = kWall/k;
//     
//     //- Fourier number
//     scalar Fo = 2.9432752*thermalDiffusivity/(radius*U_dot_nw);
//     
//     const scalar e0 = 0.87093*(p.Tsolid()-Twall)*Ac*sqrt(tc)/(pow((rhoSolid*CpSolid*k),-0.5)+pow((materialDensityWall*CpWall*kWall),-0.5));
//     
//     scalar gFo = 0.0;
//     
//     scalar Lamda = 0.0;
//     
//     if(alphaPrime > 1)
//     {
// //         if(Fo >= 0.5*sqr(radius/Rc))
// //         {
// //             //- the case of Fo -> infinity
// //             Lamda = kPrime/(1.0+kPrime)+sqrt(alphaPrime)/(1.0+kPrime);
// //         }
// //         else
// //         {
//             scalar beta = 0.95*Fo/(0.95*Fo+0.5*sqrt(Fo)+0.12);
//             
//             scalar gamma = (1.0-exp(-10.0*sqr(log(alphaPrime))))*0.0777/(0.7666*pow(Fo,-1.0/3.0)+sqrt(Fo));
//             
//             Lamda = 1.0+beta*(sqrt(alphaPrime)-1.0)/(1.0+kPrime)+gamma*4.0*kPrime/sqr(1+kPrime);
// //         }
//         
//         gFo = 0.605039+pow((0.155994+(1.182618-0.332298*exp(-2.0*Fo))*Fo),0.5);
//         
//     }
//     else if(alphaPrime < 1)
//     {
//         Fo = Fo*alphaPrime;
//         
//         alphaPrime = 1.0/alphaPrime;
//         
//         kPrime = 1.0/kPrime;
//         
//         scalar beta = 0.95*Fo/(0.95*Fo+0.5*sqrt(Fo)+0.12);
//             
//         scalar gamma = (1.0-exp(-10*sqr(log(alphaPrime))))*0.0777/(0.7666*pow(Fo,-1.0/3.0)+sqrt(Fo));
//         
//         Lamda = 1.0+beta*(sqrt(alphaPrime)-1.0)/(1.0+kPrime)+gamma*4.0*kPrime/sqr(1+kPrime);
//         
//         gFo = 0.605039+pow((0.155994+(1.182618-0.332298*exp(-2.0*alphaPrime*Fo))*alphaPrime*Fo),0.5);
//     }
//     else if(alphaPrime == 1)
//     {
//         Lamda = 1.0;
//         
//         gFo = 0.605039+pow((0.155994+(1.182618-0.332298*exp(-2.0*Fo))*Fo),0.5);
//     }
//     
//     //- corrected conductive heat transfer per impact, unit: J
//     //scalar e = Lamda*gFo*e0;
//     
//     return Lamda*gFo*e0;
    return 0.0;
    
}


Foam::scalar Foam::solidParticleCouplingCloud::particleLiquidDensityCorrection
(
    const word materialName,
    const scalar temperature
) const
{
    scalar liquidPhaseDensity = 0.0;
//     if(materialName == "Al")
//     {
//         //- Aluminum
//             
//         //- 148.3*16.0184634*(1.0 - 126e-6*9/5*(T-1679*5/9))
//     
//         liquidPhaseDensity = 2735.538122*(1.0-2.268e-4*(temperature - 932.77778));
//     }
    if(materialName == "Al2O3")
    {
        //- Aluminum Oxide
        //- 188.1*16.0184634*(1.0- 207.3e-6*9/5*(T-4188.6*5/9))
        
        liquidPhaseDensity = 3013.072966*(1.0-3.7314e-4*(temperature - 2327.0));
    }
    else if(materialName == "ZrO2")
    {
        //- the equation reference from 
        //-"Thermophysical properties of liquid UO2,ZrO2 and corium by molecular dynamics and predictive models", W.K.Kim, J.H.Shim, M. Kaviany, Journal of Nuclear Materials, Page.130, Equation (24)
        liquidPhaseDensity = (8.62 -0.89e-3*temperature)*1000;
    }
    else
    {
        FatalErrorInFunction
        << "    Unclear Material Name ! "
        << "    This material is not defined in particleLiquidDensityCorrection() !"
        << "    Please check!"
        << nl
        << exit(FatalError); 
    }
    
    return liquidPhaseDensity;
}


Foam::scalar Foam::solidParticleCouplingCloud::particleLatentHeatOfFusion
(
    const word materialName
) const
{
    scalar hf = 0.0;
    if(materialName == "Al2O3")
    {
        
        hf = 1.07e6;//J/kg
    }
    else if(materialName == "ZrO2")
    {
        //- the value reference from 
        //-"Thermophysical properties of liquid UO2,ZrO2 and corium by molecular dynamics and predictive models", W.K.Kim, J.H.Shim, M. Kaviany, Journal of Nuclear Materials, Page.132 
        //- The heat of fusion in this paper is enthalpy of fusion, which is known as latent heat of fusion
        hf = 0.26e5;//J/kg
    }
    else
    {
        FatalErrorInFunction
        << "    Unclear Material Name ! "
        << "    This material is not defined in phaseChangeModelConst() !"
        << "    Please check!"
        << nl
        << exit(FatalError); 
    }
    
    return hf;
}

Foam::scalar Foam::solidParticleCouplingCloud::phaseChangeModelConst
(
     const word materialName
) const
{
    scalar Aconst = 0.0;

    if(materialName == "Al2O3")
    {
        // - 2.7e-6 ms^-1*K^-1.8
        Aconst = 2.7e-6;//m*s^-1*K^-1.8
    }
    else if(materialName == "ZrO2")
    {
        Aconst = 0.72e-6;
    }
    else
    {
        FatalErrorInFunction
        << "    Unclear Material Name ! "
        << "    This material is not defined in phaseChangeModelConst() !"
        << "    Please check!"
        << nl
        << exit(FatalError); 
    }
    
    return Aconst;
}

// Foam::scalar Foam::solidParticleCouplingCloud::dParcel
// (
//     solidParticleCoupling& p
// )
// {
//     scalar parcelDiameter = 0;
//     if(dsmcCloudReference_->axisymmetric())
//     {
//         parcelDiameter = p.D()*pow((p.RWF()*nSolidParticles_),1/3);
//     }
//     else
//     {
//         parcelDiameter = p.D()*pow(nSolidParticles_,1/3);
//     }
//
//
//     return parcelDiameter;
// }


void Foam::solidParticleCouplingCloud::resetUcorrect()
{
    forAllIter(solidParticleCouplingCloud, *this, iter)
    {
        solidParticleCoupling* p = &iter();
        
        p->UCorrect() = Zero;
    }
}

void Foam::solidParticleCouplingCloud::updateParticleNumberSequence()
{
    //- relabel particles
    label i = 0;
    forAllIter(solidParticleCouplingCloud,*this,iter)
    {
        solidParticleCoupling& p = iter();
        p.numSeq() = i;
        i++;
    }
}

void Foam::solidParticleCouplingCloud::MPPICprocedures
(
    solidParticleCoupling::trackingData& tdSolid
)
{
    
    if(enableMPPICMethod_)
    {
        Info << "    Conducting MP-PIC method"<< endl;
        
        updateMppicCellList();
        
        //- update particle RWF 
        if(dsmcCloudReference_->axisymmetric())
        {
            axisymmetricWeighting();
            buildCellOccupancy();
        }
        
        solidParticleCoupling::TrackingData<solidParticleCouplingCloud > tdmppic(*this);
        // Damping
        if(mppicDampingModel_->active())
        {
            tdmppic.updateMPPICAverages(*this);
            mppicDampingModel_->cacheFields(true);
//             tdSolid.part() = solidParticleCoupling::trackingData::DampingNoTrack;
//             Cloud<solidParticleCoupling>::move(*this, tdSolid, mesh_.time().deltaTValue());//- just update Ucorrect
            mppicDampingModel_->calUcorrect();
            
            tdSolid.part() = solidParticleCoupling::trackingData::CorrectTrack;
            Cloud<solidParticleCoupling>::move(*this, tdSolid, mesh_.time().deltaTValue());
            mppicDampingModel_->cacheFields(false);
//             buildCellOccupancy();
            resetUcorrect();
        }
        
        // Packing
        if(mppicPackingModel_->active())
        {
            tdmppic.updateMPPICAverages(*this);
            mppicPackingModel_->cacheFields(true);
//             tdSolid.part() = solidParticleCoupling::trackingData::PackingNoTrack;
//             Cloud<solidParticleCoupling>::move(*this, tdSolid, mesh_.time().deltaTValue());//- just update Ucorrect
            mppicPackingModel_->calUcorrect();
            
            tdSolid.part() = solidParticleCoupling::trackingData::CorrectTrack;
            Cloud<solidParticleCoupling>::move(*this, tdSolid, mesh_.time().deltaTValue());
            mppicPackingModel_->cacheFields(false);
//             buildCellOccupancy();
            resetUcorrect();
        }
        
        // Isotropy
        if(mppicIsotropyModel_->active())
        {
            tdmppic.updateMPPICAverages(*this);
            mppicIsotropyModel_->calculate();
        }
        
        buildCellOccupancy();
        
        clearMppicCellList();
    }
}


/*---------- Begin of Evolve ----------*/
// -- Evolve function used in rarefiedMultiphaseFoam.C
void Foam::solidParticleCouplingCloud::evolve()
{
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();//****
    
    solidParticleCoupling::trackingData tdSolid(*this);
    
    boundaries_.controlBeforeCollisions();//- empty function
    
    //- interphase calculation, do not include velocity update in it
    interphaseCouplingCal();

    //- upadte particle velocity half time-step due to aerodynamic force
    updateVelocityHalfTimeStep();
    
    //- include the influence of gravitaional field, if activated
    controllers_.controlBeforeMove();
    //- sub-functions of boundaries_ are in dsmcSolidboundaries.C
    boundaries_.controlAfterCollisions();//- empty function
    
    //- release solid particles 
    boundaries_.controlBeforeMove();
    
    // Move the particles ballistically with their average velocities
    //- position is updated (without MPPIC)
//     particleMPPICMove_ = false;
    tdSolid.part() = solidParticleCoupling::trackingData::LinearTrack;

    Cloud<solidParticleCoupling>::move(*this, tdSolid, mesh_.time().deltaTValue());

    buildCellOccupancy();
    
    if(solidCollisionDetectionModel_->active())
    {
        updateVelocityHalfTimeStepBack();
    
        controllers_.controlBeforeCollisions();//- update velocity because of grvity
    
        particleParticleCollisions();
        
        updateVelocityFullTimeStep();
        
        controllers_.controlAfterCollisions();
    }
    else
    {
        updateVelocityHalfTimeStep();
        controllers_.controlAfterMove();
    }
    
    //- MPPIC method
    MPPICprocedures(tdSolid);
    
    if(dsmcCloudReference_->axisymmetric())
    {
        axisymmetricWeighting();
        buildCellOccupancy();
    }
    
    updateParticleVolumeFraction();
    fields_.calculateFields();
    fields_.writeFields();
    
    controllers_.calculateProps();//- empty for solidGravitationalAccelerationController
    controllers_.outputResults();//- empty for solidGravitationalAccelerationController
    
    boundaries_.calculateProps();
    boundaries_.outputResults();

    solidBoundaryMeas_.clean();
    solidCellMeas_.clean();
    
    info();
    
}
/*---------- End of Evolve ----------*/

void Foam::solidParticleCouplingCloud::axisymmetricWeighting()
{
    const auto& cellOccupancy = this->cellOccupancy();

    forAll(cellOccupancy, c)
    {
        for (solidParticleCoupling* p : cellOccupancy[c])
        {
            const point& cC = mesh_.cellCentres()[c];
            const scalar radius = cC.y();
                        
            const scalar oldRadialWeight = p->RWF();
            const scalar newRadialWeight = 1.0 + dsmcCloudReference()->maxRWF()*(radius/dsmcCloudReference()->radialExtent());
//             const scalar newRadialWeight = dsmcCloudReference()->axiRWF(cC);
            
                
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

                    addNewSolidParticle
                    (
                        position,
                        Us,
                        p->omega(),
                        p->UCorrect(),
                        p->F(),
                        p->D(),
                        p->T(),
                        p->RWF(),
                        p->CzRatio(),
                        c,
                        -1,
                        p->typeID(),
                        p->phaseState(),
                        p->newSolidParticle()
                        
                    );
                        
                    prob -= 1.0;
                }
                    
                if(prob > rndGenS_.sample01<scalar>())
                {
                    const vector position = p->position();
                    
                    vector Us = p->U();
                    
                    Us.z() *= -1.0;

                    addNewSolidParticle
                    (
                        position,
                        Us,
                        p->omega(),
                        p->UCorrect(),
                        p->F(),
                        p->D(),
                        p->T(),
                        p->RWF(),
                        p->CzRatio(),
                        c,
                        -1,
                        p->typeID(),
                        p->phaseState(),
                        p->newSolidParticle()
                    );
                }
            }
            
            if(newRadialWeight > oldRadialWeight)
            {           
                //particle might be deleted
                if((oldRadialWeight/newRadialWeight) < rndGenS_.sample01<scalar>())
                {
                    deleteParticle(*p);
                } 
            } 
        }
    }
}

void Foam::solidParticleCouplingCloud::addNewSolidParticle
(
    const vector& position,
    const vector& U,
    const vector& omega,
    const vector& UCorrect,
    const vector& F,
    const scalar D,
    const scalar T,
    const scalar RWF,
    const scalar CzRatio,
    const label cellI,
    const label numSeq,
    const label typeID,
    const label phaseState,
    const label newSolidParticle
)
{
    //- turn particle diameter to parcel diameter
//     D = D*pow(nSolidParticles_,1/3);

    solidParticleCoupling* pPtr = new solidParticleCoupling
    (
        mesh_,
        position,
        U,
        omega,
        UCorrect,
        F,
        D,
        T,
        RWF,
        CzRatio,
        cellI,
        numSeq,
        typeID,
        phaseState,
        newSolidParticle
    );

    addParticle(pPtr);
}

void Foam::solidParticleCouplingCloud::info() const
{
    label nSolidSimulators = this->size();
    reduce(nSolidSimulators, sumOp<label>());

    
    Info << "    Number of solid particles       = " << nSolidSimulators<< endl;
    
    Info << endl;
    
}

// ************************************************************************* //
