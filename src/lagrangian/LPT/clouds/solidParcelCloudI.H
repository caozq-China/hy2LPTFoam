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

#include "fvmSup.H"
#include "SortableList.H"
#include "constants.H"
using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const Foam::word& Foam::solidParcelCloud::cloudName() const
{
    return cloudName_;
}

inline const Foam::label& Foam::solidParcelCloud::interiorInteractionType() const
{
    return interiorInteractionType_;
}

// inline const Foam::dynamicFvMesh& Foam::solidParcelCloud::mesh() const
inline const Foam::fvMesh& Foam::solidParcelCloud::mesh() const
{
    return mesh_;
}

inline const Foam::IOdictionary& Foam::solidParcelCloud::particleProperties() const
{
    return particleProperties_;
}

inline const Foam::cloudSolution& Foam::solidParcelCloud::solution() const
{
    return solution_;
}

inline const Foam::dictionary& Foam::solidParcelCloud::subModelProperties() const
{
    return subModelProperties_;
}


inline const Foam::List<Foam::word>& Foam::solidParcelCloud::typeIdList() const
{
    return typeIdList_;
}

inline Foam::scalar Foam::solidParcelCloud::nParticle() const
{
    // return coordSystem().nParticlesOrg();
    return nParticle_;
}

inline Foam::scalar Foam::solidParcelCloud::radialExtent() const
{
    return radialExtent_;
}

inline Foam::scalar Foam::solidParcelCloud::maxRWF() const
{
    return maxRWF_;
}

inline bool Foam::solidParcelCloud::axisymmetric() const
{
    return axisymmetric_;
}

//inline const Foam::volScalarField& Foam::solidParcelCloud::nParticles() const
//{
//    return coordSystem().nParticles();
//}


//inline Foam::scalar Foam::solidParcelCloud::nParticles(const label celli) const
//{
//    return coordSystem().nParticles(celli);
//}


//inline Foam::scalar Foam::solidParcelCloud::nParticles
//(
//    const label patchi,
//    const label facei
//) const
//{
//    return coordSystem().nParticles(patchi, facei);
//}

inline Foam::List<Foam::DynamicList<Foam::solidParcel*>>& Foam::solidParcelCloud::cellOccupancy()
{
    if(cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}

inline const Foam::scalarField& Foam::solidParcelCloud::cellLengthScale() const
{
    return cellLengthScale_;
}

inline const Foam::List<Foam::solidParcel::constantProperties>&
Foam::solidParcelCloud::constProps() const
{
    return constProps_;
}

inline const Foam::rho2ReactionThermo& Foam::solidParcelCloud::thermo() const
{
    return thermo_;
}

inline const Foam::boundaryMeasurements& Foam::solidParcelCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}

inline Foam::boundaryMeasurements& Foam::solidParcelCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}

inline const Foam::solidFieldProperties& Foam::solidParcelCloud::fieldProps() const
{
    return fieldProps_;
}

inline Foam::solidFieldProperties& Foam::solidParcelCloud::fieldProps()
{
    return fieldProps_;
}

inline const Foam::diffusionBasedSmoothing& Foam::solidParcelCloud::DBS() const
{
    return DBS_;
}

inline Foam::diffusionBasedSmoothing& Foam::solidParcelCloud::DBS()
{
    return DBS_;
}

inline const Foam::diffusionBasedGasFiltering& Foam::solidParcelCloud::DBGF() const
{
    return DBGF_;
}

inline Foam::diffusionBasedGasFiltering& Foam::solidParcelCloud::DBGF()
{
    return DBGF_;
}

inline const Foam::PackingModel& Foam::solidParcelCloud::packingModel() const
{
    return packingModel_;
}

inline Foam::PackingModel& Foam::solidParcelCloud::packingModel()
{
    return packingModel_();
}

inline const Foam::DampingModel& Foam::solidParcelCloud::dampingModel() const
{
    return dampingModel_;
}

inline Foam::DampingModel& Foam::solidParcelCloud::dampingModel()
{
    return dampingModel_();
}

inline const Foam::IsotropyModel& Foam::solidParcelCloud::isotropyModel() const
{
    return isotropyModel_;
}

inline Foam::IsotropyModel& Foam::solidParcelCloud::isotropyModel()
{
    return isotropyModel_();
}

inline const Foam::HeatTransferModel& Foam::solidParcelCloud::heatTransfer() const
{
    return heatTransferModel_;
}

inline Foam::HeatTransferModel& Foam::solidParcelCloud::heatTransfer()
{
    return heatTransferModel_();
}

inline const Foam::solidParcel::constantProperties& Foam::solidParcelCloud::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorIn("Foam::solidParcelCloud::constProps(label typeId)")
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }
    return constProps_[typeId];
}

inline Foam::Random& Foam::solidParcelCloud::rndGen()
{
    return rndGen_;
}

inline Foam::volScalarField& Foam::solidParcelCloud::theta()
{
    return theta_;
}

inline const Foam::volVectorField& Foam::solidParcelCloud::UFilter() const
{
    return UFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::rhoFilter() const
{
    return rhoFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::muFilter() const
{
    return muFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::TtrFilter() const
{
    return TtrFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::CpFilter() const
{
    return CpFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::kappaFilter() const
{
    return kappaFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::MaFilter() const
{
    return MaFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::gammaFilter() const
{
    return gammaFilter_;
}

inline const Foam::volScalarField& Foam::solidParcelCloud::omegaFilter() const
{
    return omegaFilter_;
}

inline const Foam::dimensionedVector& Foam::solidParcelCloud::g() const
{
    return g_;
}

inline const Foam::forceList& Foam::solidParcelCloud::forces() const
{
    return forcesList_;
}

inline Foam::forceList& Foam::solidParcelCloud::forces()
{
    return forcesList_;
}

inline const Foam::solidBoundaries& Foam::solidParcelCloud::boundaries() const
{
    return boundaries_;
}

inline Foam::solidBoundaries& Foam::solidParcelCloud::boundaries()
{
    return boundaries_;
}

//inline const Foam::coordinateSystemType& Foam::solidParcelCloud::coordSystem() const
//{
//    return particleCoordinateSystem_;
//}

//inline Foam::coordinateSystemType& Foam::solidParcelCloud::coordSystem()
//{
//    return particleCoordinateSystem_();
//}

inline const Foam::integrationScheme& Foam::solidParcelCloud::UIntegrator() const
{
    return UIntegrator_();
}

inline const Foam::integrationScheme& Foam::solidParcelCloud::TIntegrator() const
{
    return TIntegrator_();
}

inline bool Foam::solidParcelCloud::PreduceD() const
{
    return PreduceD_;
}


inline Foam::DimensionedField<Foam::vector, Foam::volMesh>& Foam::solidParcelCloud::momentumTrans()
{
    return momentumTrans_();
}

inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>& Foam::solidParcelCloud::momentumTrans() const
{
    return momentumTrans_();
}

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::momentumCoeff()
{
    return momentumCoeff_();
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::momentumCoeff() const
{
    return momentumCoeff_();
}

inline Foam::tmp<Foam::fvVectorMatrix> Foam::solidParcelCloud::SU(volVectorField& U) const
{
    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const volScalarField::Internal
                Vdt(mesh_.V()*this->db().time().deltaT());

            return 
                momentumTrans()/Vdt - fvm::Sp(momentumCoeff()/Vdt, U) + momentumCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm.ref();

            fvm.source() = -momentumTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::qTrans()
{
    return qTrans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::qTrans() const
{
    return qTrans_();
}

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::qCoeff()
{
    return qCoeff_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::qCoeff() const
{
    return qCoeff_();
}

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::wTrans()
{
    return wTrans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>& Foam::solidParcelCloud::wTrans() const
{
    return wTrans_();
}

inline Foam::tmp<Foam::fvScalarMatrix> Foam::solidParcelCloud::SE(volScalarField& e) const
{

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("e"))
        {
            const volScalarField Cv(thermo_.Cv_t());

            const volScalarField::Internal
                Vdt(this->mesh().V()*this->db().time().deltaT());

            return
                
            qTrans()/Vdt
              - fvm::SuSp(qCoeff()/(Cv*Vdt), e)
              + qCoeff()/(Cv*Vdt)*e;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(e, dimEnergy/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -( qTrans()/(this->db().time().deltaT()) + wTrans() );

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(e, dimEnergy/dimTime));
}

//- check

inline Foam::scalar Foam::solidParcelCloud::Tmax() const
{
    scalar T = -GREAT;
    scalar n = 0;
    forAllConstIter(solidParcelCloud, *this, iter)
    {
        const solidParcel& p = iter();
        T = max(T, p.T());
        n++;
    }

    reduce(T, maxOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}

inline Foam::scalar Foam::solidParcelCloud::Tmin() const
{
    scalar T = GREAT;
    scalar n = 0;
    forAllConstIter(solidParcelCloud, *this, iter)
    {
        const solidParcel& p = iter();
        T = min(T, p.T());
        n++;
    }

    reduce(T, minOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


inline Foam::label Foam::solidParcelCloud::nParcelsInSystem() const
{
    return this->size();
}

inline const Foam::tmp<Foam::volScalarField> Foam::solidParcelCloud::np() const
{
    tmp<volScalarField> tnp
    (
        new volScalarField
        (
            IOobject
            (
                "nRealParticles",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& np = tnp.ref();
    forAllConstIter(solidParcelCloud, *this, iter)
    {
        const solidParcel& p = iter();
        const label celli = p.cell();

        np[celli] += nParticle()*p.RWF();
    }

    return tnp;
}

inline void Foam::solidParcelCloud::updateTheta()
{
    theta_.ref().field() = 0.0;

    forAllConstIter(solidParcelCloud, *this, iter)
    {
        const solidParcel& p = iter();
        const label celli = p.cell();

        theta_.ref()[celli] += nParticle()*p.RWF()*p.volume();
    }

    theta_.primitiveFieldRef() /= mesh_.V();
    theta_.correctBoundaryConditions();
}


inline Foam::scalar Foam::solidParcelCloud::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(solidParcelCloud, *this, iter)
    {
         const solidParcel& p = iter();

        sysMass += nParticle()*p.RWF()*p.mass();
    }

    return sysMass;
}


inline Foam::vector Foam::solidParcelCloud::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(solidParcelCloud, *this, iter)
    {
        const solidParcel& p = iter();

        linearMomentum += nParticle()*p.RWF()*p.mass()*p.U();
    }

    return linearMomentum;
}


inline Foam::scalar Foam::solidParcelCloud::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(solidParcelCloud, *this, iter)
    {
        const solidParcel& p = iter();

        linearKineticEnergy += nParticle()*p.RWF()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}

// ************************************************************************* //
