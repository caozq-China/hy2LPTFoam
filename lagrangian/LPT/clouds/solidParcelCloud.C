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

#include "solidParcelCloud.H"
#include "zeroGradientFvPatchFields.H"
#include "constants.H"
#include "interpolation.H"
#include "subCycleTime.H"
#include "HeatTransferModel.H"
using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<solidParcel>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::solidParcelCloud::setModels()
{
    if(interiorInteractionType_ == 0)
    {
        Info<<"Attention: No particle-particle interaction model is activated!"<<endl;
    }
    else if(interiorInteractionType_ == 1)
    {
    }
    else if(interiorInteractionType_ == 2)
    {
        Info<<"Attention:  Multiphase particle-in-cell(MPPIC) method is activated!"<<endl;
        packingModel_ = PackingModel::New(subModelProperties_.subDict("MPPIC"),*this);
        dampingModel_ = DampingModel::New(subModelProperties_.subDict("MPPIC"),*this);
        isotropyModel_= IsotropyModel::New(subModelProperties_.subDict("MPPIC"),*this);

    }
    else
    {
        FatalErrorInFunction
            << "    Warning! "
            << "    Unknown particle-particle interaction model is used!"
            << "    If no particle-particle interaction is considered, "
            << "    please put 0 in interiorInteractionType"
            << nl
            << exit(FatalError);
    }

    UIntegrator_.reset
    (
        integrationScheme::New
        (
            "U",
            solution_.integrationSchemes()
        ).ptr()
    );

    TIntegrator_.reset
    (
        integrationScheme::New
        (
            "T",
            solution_.integrationSchemes()
        ).ptr()
    );

    

    heatTransferModel_ = HeatTransferModel::New
    (
        subModelProperties_,
        *this
    );
}

void Foam::solidParcelCloud::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary solidProperties
    (
        particleProperties_.subDict("constantProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& particleDict(solidProperties.subDict(id));

        constProps_[i] = solidParcel::constantProperties(particleDict);
    }
}

void Foam::solidParcelCloud::buildCellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        cellOccupancyPtr_.reset
        (
            new List<DynamicList<solidParcel*>>(mesh_.nCells())
        );
    }
    else if (cellOccupancyPtr_().size() != mesh_.nCells())
    {
        // If the size of the mesh has changed, reset the
        // cellOccupancy size
        cellOccupancyPtr_().setSize(mesh_.nCells());
    }

    List<DynamicList<solidParcel*>>& cellOccupancy = cellOccupancyPtr_();

    forAll(cellOccupancy, celli)
    {
        cellOccupancy[celli].clear();
    }

    forAllIter(solidParcelCloud, *this, iter)
    {
        cellOccupancy[iter().cell()].append(&iter());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//-Constructer of running simulations
Foam::solidParcelCloud::solidParcelCloud
(
    Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& Ma,
    const volScalarField& gamma,
    const volScalarField& omega,
    const dimensionedVector& g,
    const rho2ReactionThermo& thermo,
    bool readFields
)
:
    Cloud<solidParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName+"Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    subModelProperties_(particleProperties_.subDict("subModels")),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    constProps_(),
    g_(g),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    radialExtent_(0.0),
    maxRWF_(1.0),
    interiorInteractionType_(readLabel(particleProperties_.lookup("interiorInteractionType"))),
    axisymmetric_(Switch(particleProperties_.lookup("axisymmetric"))),
    PreduceD_(particleProperties_.lookupOrDefault("PreduceD",false)),
    hasWallImpactDistance_(particleProperties_.lookupOrDefault("hasWallImpactDistance",false)),
    thermo_(thermo),
    U_(U),
    Ma_(Ma),
    gamma_(gamma),
    omega_(omega),
    UFilter_
    (
        IOobject
        (
            "UFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_
    ),
    rhoFilter_
    (
        IOobject
        (
            "rhoFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.rho()
    ),
    muFilter_
    (
        IOobject
        (
            "muFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.mu()
    ),
    TtrFilter_
    (
        IOobject
        (
            "TtrFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.T()
    ),
    CpFilter_
    (
        IOobject
        (
            "CpFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.Cp_t()
    ),
    kappaFilter_
    (
        IOobject
        (
            "kappaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.kappa()
    ),
    MaFilter_
    (
        IOobject
        (
            "MaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ma_
    ),
    gammaFilter_
    (
        IOobject
        (
            "gammaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        gamma_
    ),
    omegaFilter_
    (
        IOobject
        (
            "omegaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        omega_
    ),
    theta_
    (
        IOobject
        (
            "theta",
            t.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    np_
    (
        IOobject
        (
            "nrp",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    cellLengthScale_(cbrt(mesh_.V())),
    //particleCoordinateSystem_(coordinateSystemType::New(t,mesh,*this)),
    cellOccupancyPtr_(),
    rndGen_(label(clock::getTime()) + 7183*Pstream::myProcNo()),
    solution_(mesh_, particleProperties_.subDict("solution")),
    forcesList_
    (
        *this,
        mesh_,
        subModelProperties_.subDict("Forces"),
        true
    ),
    boundaries_(t,mesh,*this),
    boundaryMeas_(mesh, *this, true),
    fieldProps_(t, mesh, *this),
    DBS_(particleProperties_, t, mesh_, *this),
    DBGF_(particleProperties_, t, mesh_, *this),
    UIntegrator_(nullptr),
    TIntegrator_(nullptr),
    packingModel_(),
    dampingModel_(),
    isotropyModel_(),
    heatTransferModel_(nullptr),
    momentumTrans_
    (
        new volVectorField::Internal
        // new volVectorField
        (
            IOobject
            (
                "momentumTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimMass*dimVelocity, Zero)
        )
    ),
    momentumCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "momentumCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimMass, 0.0)
        )
    ),
    qTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "qTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy, 0.0)
        )
    ),
    qCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "qCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTemperature, 0.0)
        )
    ),
    wTrans_
    (
        new volScalarField::Internal
        // new volScalarField
        (
            IOobject
            (
                "wTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
        )
    )
{   
    setModels();

    if (readFields)
    {
        solidParcel::readFields(*this);
    }
    
    resetSourceTerms();

    buildConstProps();
    if(axisymmetric_)
    {
        radialExtent_ = 
                readScalar(particleProperties_.lookup("radialExtentOfDomain"));
        maxRWF_ = 
            readScalar(particleProperties_.lookup("maxRadialWeightingFactor"));
        maxRWF_ -= 1.0;
    }
    //coordSystem().checkCoordinateSystemInputs();

    buildCellOccupancy();
    
    fieldProps_.createFields();
    boundaries_.setInitialConfig();
}

//- Constructer of Initialization
Foam::solidParcelCloud::solidParcelCloud
(
    Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& Ma,
    const volScalarField& gamma,
    const volScalarField& omega,
    const dimensionedVector& g,
    const rho2ReactionThermo& thermo,
    const IOdictionary& solidInitialiseDict,
    const bool& clearFields
)
:
    Cloud<solidParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName+"Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    subModelProperties_(particleProperties_.subDict("subModels")),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    constProps_(),
    g_(g),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    radialExtent_(0.0),
    maxRWF_(1.0),
    interiorInteractionType_(readLabel(particleProperties_.lookup("interiorInteractionType"))),
    axisymmetric_(Switch(particleProperties_.lookup("axisymmetric"))),
    PreduceD_(false),
    hasWallImpactDistance_(false),
    thermo_(thermo),
    U_(U),
    Ma_(Ma),
    gamma_(gamma),
    omega_(omega),
    UFilter_
    (
        IOobject
        (
            "UFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_
    ),
    rhoFilter_
    (
        IOobject
        (
            "rhoFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.rho()
    ),
    muFilter_
    (
        IOobject
        (
            "muFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.mu()
    ),
    TtrFilter_
    (
        IOobject
        (
            "TtrFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.T()
    ),
    CpFilter_
    (
        IOobject
        (
            "CpFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.Cp_t()
    ),
    kappaFilter_
    (
        IOobject
        (
            "kappaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.kappa()
    ),
    MaFilter_
    (
        IOobject
        (
            "MaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ma_
    ),
    gammaFilter_
    (
        IOobject
        (
            "gammaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        gamma_
    ),
    omegaFilter_
    (
        IOobject
        (
            "omegaFilter",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        omega_
    ),
    theta_
    (
        IOobject
        (
            "theta",
            t.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    np_
    (
        IOobject
        (
            "nrp",
            t.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    cellLengthScale_(cbrt(mesh_.V())),
    //particleCoordinateSystem_(coordinateSystemType::New(t,mesh,*this)),
    cellOccupancyPtr_(),
    rndGen_(label(clock::getTime()) + 7183*Pstream::myProcNo()),
    solution_(mesh_, particleProperties_.subDict("solution")),
    forcesList_
    (
        *this,
        mesh_,
        subModelProperties_.subDict("Forces"),
        true
    ),
    boundaries_(t,mesh,*this),
    boundaryMeas_(mesh, *this),
    fieldProps_(t, mesh),
    DBS_(particleProperties_, t, mesh_, *this),
    DBGF_(particleProperties_, t, mesh_, *this),
    UIntegrator_(nullptr),
    TIntegrator_(nullptr),
    packingModel_(),
    dampingModel_(),
    isotropyModel_(),
    heatTransferModel_(nullptr),
    momentumTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                "momentumTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimMass*dimVelocity, Zero)
        )
    ),
    momentumCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "momentumCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimMass, 0.0)
        )
    ),
    qTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "qTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy, 0.0)
        )
    ),
    qCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "qCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTemperature, 0.0)
        )
    ),
    wTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "wTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
        )
    )
{
    if(!clearFields)
    {
        solidParcel::readFields(*this);
    }

    label initialSolidParticles = this->size();

    if (Pstream::parRun())
    {
        reduce(initialSolidParticles, sumOp<label>());
    }

    if(clearFields)
    {
        Info << "Clearing existing field of solid particles " << endl;

        IDLList<solidParcel>::clear();

        initialSolidParticles = 0;
    }

    if(axisymmetric_)
    {
        radialExtent_ = 
            readScalar(particleProperties_.lookup("radialExtentOfDomain"));
        maxRWF_ = 
            readScalar(particleProperties_.lookup("maxRadialWeightingFactor"));
        maxRWF_ -= 1.0;
    }
    
    buildConstProps();

    //- --Begin initialization-- -//
    allConfigurations solidIntialization(solidInitialiseDict, *this);
    solidIntialization.setInitialConfig();//- pointer pointing to setInitialConfiguration(), e.g. zoneFill  
    //- --End initialization-- -//

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
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidParcelCloud::~solidParcelCloud()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::labelList Foam::solidParcelCloud::getTypeIDs(const dictionary& dict) const
{
    const wordHashSet names(dict.lookup("typeIds"));//- [P1,P2]

    if (names.empty())
    {
        FatalErrorInFunction
            << "Entry solid phase typeIds cannot be empty in " << dict.name()
            << exit(FatalError);
    }

    labelList typeIDs(names.size(), -1);//- [-1,-1]

    label i = 0;
    forAllIters(names, iter)
    {
        const word& name = iter();

        label id = findIndex(typeIdList_, name);

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
// = = = = = = = = = = = = = = = = = = = = =
void Foam::solidParcelCloud::evolveCloud
(
    solidParcel::trackingData& td
)
{
    boundaries_.controlBeforeCollisions();//- empty function

    //- general boundary will release solid particles
    boundaries_.controlBeforeMove();

    td.part() = solidParcel::trackingData::LinearTrack;
    Cloud<solidParcel>::move(td,mesh_.time().deltaTValue());//- move particle
    if(axisymmetric_)
    {
        axisymmetricWeighting();
        buildCellOccupancy();
    }
    else
    {
        buildCellOccupancy();// update cell occupancy
    }
    //- Radial weighting for non-Cartesian flows (e.g., axisymmetric). This is
    // where parcels will receive their new RWF and will possibly be cloned or
    // deleted.
    //coordSystem().evolve();

    if(interiorInteractionType_==1)
    {   

    }
    else if(interiorInteractionType_==2)//- MPPIC
    {
        td.cloud().MPPICMotion(td);
    }

    boundaries_.controlAfterCollisions();//- for general Boundary condition

    // update Volume Fraction
    updateTheta();
}

void Foam::solidParcelCloud::evolve()
{   
    label nGeometricD = mesh_.nGeometricD();

    Info<< "\nSolving "<< nGeometricD << "-D"<<" lagrangian cloud " << endl;

    solidParcel::trackingData td(*this,solidParcel::trackingData::LinearTrack);

    forcesList_.cacheFields(true);

    buildCellOccupancy();

    resetSourceTerms();

    //- end of preEvolve()

    evolveCloud(td);

    fieldProps_.calculateFields();
    fieldProps_.writeFields();

    forcesList_.cacheFields(false);
    boundaryMeas_.clean();
    
    td.cloud().info();
}

void Foam::solidParcelCloud::MPPICMotion(solidParcel::trackingData& td)
{
    forces().setCalcNonCoupled(false);
    forces().setCalcCoupled(false);

    if(dampingModel_->active() || packingModel_->active())
    {
        td.updateAverages(*this);

        if(dampingModel_->active())
        {
            dampingModel_->cacheFields(true);
            td.part() = solidParcel::trackingData::DampingNoTrack;
            Cloud<solidParcel>::move(td,mesh_.time().deltaTValue());//- move
            td.part() = solidParcel::trackingData::CorrectTrack;
            Cloud<solidParcel>::move(td,mesh_.time().deltaTValue());//- move
        }
        if(packingModel_->active())
        {
            packingModel_->cacheFields(true);
            td.part() = solidParcel::trackingData::PackingNoTrack;
            Cloud<solidParcel>::move(td,mesh_.time().deltaTValue());//- move
            td.part() = solidParcel::trackingData::CorrectTrack;
            Cloud<solidParcel>::move(td,mesh_.time().deltaTValue());//- move
        }
    }

    // td.part() = solidParcel::trackingData::LinearTrack;
    // Cloud<solidParcel>::move(td,mesh_.time().deltaTValue());//- move particle

    if(isotropyModel_->active())
    {
        td.updateAverages(*this);
        isotropyModel_->calculate();
    }

    //buildCellOccupancy();// update cell occupancy
    if(axisymmetric_)
    {
        axisymmetricWeighting();
        buildCellOccupancy();
    }
    else
    {
        buildCellOccupancy();// update cell occupancy
    }

    //- Radial weighting for non-Cartesian flows (e.g., axisymmetric). This is
    // where parcels will receive their new RWF and will possibly be cloned or
    // deleted.
    //coordSystem().evolve();

    // forces().setCalc(true);
    forces().setCalcNonCoupled(true);
    forces().setCalcCoupled(solution().coupled());

    if(dampingModel_->active())
    {
        dampingModel_->cacheFields(false);
    }
    if(packingModel_->active())
    {
        packingModel_->cacheFields(false);
    }
}

bool Foam::solidParcelCloud::hasWallImpactDistance() const
{
    // return true;
    return hasWallImpactDistance_;
}

void Foam::solidParcelCloud::setParcelThermoProperties
(
    solidParcel& parcel,
    const scalar lagrangianDt
)
{
    parcel.rho() = constProps(parcel.typeId()).rho0();
    parcel.T() = constProps(parcel.typeId()).T0();
    parcel.Cp() = constProps(parcel.typeId()).Cp0();
}

void Foam::solidParcelCloud::updateGasProperties()
{
    UFilter_ = U_;
    rhoFilter_ = thermo_.rho();
    TtrFilter_ = thermo_.T();
    muFilter_ = thermo_.mu();
    kappaFilter_ = thermo_.kappa();
    CpFilter_ = thermo_.Cp();
    MaFilter_ = Ma_;
    gammaFilter_ = gamma_;
    omegaFilter_ = omega_;
}

void Foam::solidParcelCloud::resetSourceTerms()
{
    momentumTrans().field() = Zero;
    momentumCoeff().field() = 0.0;

    qTrans_->field() = 0.0;
    qCoeff_->field() = 0.0;

    wTrans_->field() = 0.0;

}

void Foam::solidParcelCloud::smoothSourceTerms()
{
    momentumTrans().field() *= (1-theta_.primitiveFieldRef());
    qTrans().field() *= (1-theta_.primitiveFieldRef());
    wTrans().field() *= (1-theta_.primitiveFieldRef());

    DBS().diffusion(momentumTrans_);
    DBS().diffusion(qTrans_);
    DBS().diffusion(wTrans_);

    momentumTrans().field() /= (1- theta_.primitiveFieldRef());
    qTrans().field() /= (1- theta_.primitiveFieldRef());
    wTrans().field() /= (1- theta_.primitiveFieldRef());

    if(solution().semiImplicit("e"))
    {
        qCoeff().field() *= (1-theta_.primitiveFieldRef());
        DBS().diffusion(qCoeff_);
        qCoeff().field() /= (1- theta_.primitiveFieldRef());
    }

    if(solution().semiImplicit("U"))
    {
        momentumCoeff().field() *= (1-theta_.primitiveFieldRef());
        DBS().diffusion(momentumCoeff_);
        momentumCoeff().field() /= (1- theta_.primitiveFieldRef());
    }
}

void Foam::solidParcelCloud::updateCellOccupancy()
{
    buildCellOccupancy();
}

void Foam::solidParcelCloud::filteringGasProperties()
{
    if (!solution().useGasFilterImplicit())
    {
        np_ = np();
        np_.correctBoundaryConditions();
        DBGF_.preExplicit(np_);
    }
    
    UFilter_ = DBGF_.filteredField(U_);
    rhoFilter_ = DBGF_.filteredField(thermo_.rho());
    TtrFilter_ = DBGF_.filteredField(thermo_.T());
    muFilter_ = DBGF_.filteredField(thermo_.mu());
    kappaFilter_ = DBGF_.filteredField(thermo_.kappa());
    CpFilter_ = DBGF_.filteredField(thermo_.Cp());
    MaFilter_ = DBGF_.filteredField(Ma_);
    gammaFilter_ = DBGF_.filteredField(gamma_);
    omegaFilter_ = DBGF_.filteredField(omega_);
}


void Foam::solidParcelCloud::addNewParcel
(
    const polyMesh& owner,
    const solidParcel::constantProperties& constProps,
    const vector& position,
    const vector& U0,
    const scalar RWF,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const label newParcel
)
{
    //- turn particle diameter to parcel diameter

    solidParcel* pPtr = new solidParcel
    (
        mesh_,
        constProps,
        position,
        U0,
        RWF,
        celli,
        tetFacei,
        tetPti,
        typeId,
        newParcel
    );

    addParticle(pPtr);
}

void Foam::solidParcelCloud::addNewParcel
(
    const polyMesh& owner,
    const solidParcel::constantProperties& constProps,
    const vector& position,
    const vector& U0,
    const scalar T,
    const scalar RWF,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const label newParcel
)
{
    //- turn particle diameter to parcel diameter

    solidParcel* pPtr = new solidParcel
    (
        mesh_,
        constProps,
        position,
        U0,
        T,
        RWF,
        celli,
        tetFacei,
        tetPti,
        typeId,
        newParcel
    );

    addParticle(pPtr);
}

void Foam::solidParcelCloud::axisymmetricWeighting()
{
    List<DynamicList<solidParcel*>>& cellOccupancy = cellOccupancyPtr_();
    forAll(cellOccupancy, c)
    {
        const DynamicList<solidParcel*>& pInCell = cellOccupancy[c];

        point cC = mesh_.cellCentres()[c];
        scalar radius = cC.y();

        forAll(pInCell, pIC)
        {
            solidParcel* p = pInCell[pIC];
            scalar oldRWF = p->RWF();
            scalar newRWF = 1.0;
            newRWF = 1.0 + maxRWF_*(radius/radialExtent_);
            p->RWF() = newRWF;

            if(oldRWF > newRWF)
            {
                scalar prob = (oldRWF/newRWF) - 1.0;
                while(prob>1.0)
                {
                    //add a particle and reduce prob by 1.0
                   
                    vector position = p->position();
                    
                    label cell = -1;
                    label tetFace = -1;
                    label tetPt = -1;

                    mesh_.findCellFacePt
                    (
                        position,
                        cell,
                        tetFace,
                        tetPt
                    );

                    addNewParcel
                    (
                        mesh_,
                        constProps(p->typeId()),
                        position,
                        p->U(),
                        p->RWF(),
                        cell,
                        tetFace,
                        tetPt,
                        p->typeId(),
                        p->newParcel()
                    );
                    
                    prob -= 1.0;
                }

                if(prob > rndGen_.sample01<scalar>())
                {
                    vector position = p->position();
                    
                    label cell = -1;
                    label tetFace = -1;
                    label tetPt = -1;

                    mesh_.findCellFacePt
                    (
                        position,
                        cell,
                        tetFace,
                        tetPt
                    );

                    addNewParcel
                    (
                        mesh_,
                        constProps(p->typeId()),
                        position,
                        p->U(),
                        p->RWF(),
                        cell,
                        tetFace,
                        tetPt,
                        p->typeId(),
                        p->newParcel()
                    );
                }
            }
            if(newRWF > oldRWF)
            {
                if((oldRWF/newRWF) < rndGen_.sample01<scalar>())
                {
                    deleteParticle(*p);
                }
            }
        }
    }
}

void Foam::solidParcelCloud::info()
{
    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    Info<< "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << "    Linear kinetic energy           = "
        << linearKineticEnergy << nl
        // << "    UPF                             = "
        // << UPF << nl
        // << "    DPF                             = "
        // << DPF << nl
        // << "    Axial Mass of Center            = "
        // << axialCenterOfMassOfSystem()/massInSystem() << nl
        <<endl;
}

// ************************************************************************* //
