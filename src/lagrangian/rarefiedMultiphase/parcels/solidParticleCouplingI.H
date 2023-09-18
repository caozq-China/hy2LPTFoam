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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
inline Foam::solidParticleCoupling::constantProperties::constantProperties()
:
    d_(0),
    rho_(0)
{}


inline Foam::solidParticleCoupling::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    d_(dict.get<scalar>("Diameter")),
    epsilon_(dict.get<scalar>("epsilon")),
    alpha_(dict.get<scalar>("alpha")),
    rho_(dict.get<scalar>("rho")),
    Cp_(dict.get<scalar>("specificHeat")),
    mu_(dict.get<scalar>("mu")),
    tau_(dict.get<scalar>("tau")),
    k_(0.0),
    poissonRatio_(0.0),
    elasticModuli_(0.0),
    Tm_(0.0),
    Tf_(0.0),
    effD_(0.0),
    effRho_(0.0),
    nonSperical_(dict.get<Switch>("nonSphericalParticle"))
{
    if(nonSperical_ == false)
    {
        //- if solid particles are spherical, then do nothing
        effD_ =  d_;
        effRho_ = rho_;
    }
    else if(nonSphericalSuperficalArea_ > ROOTVSMALL && nonSperical_ == true)
    {
    	Info<<"Enable nonspherical model!"<<endl;
    	
    	const dictionary& nonSphericalModelDict
        (
            dict.subDict("nonSphericalModelProPerties")
        );
        
        nonSphericalVolume_ = nonSphericalModelDict.get<scalar>("nonsphericalParticleVolume");
        nonSphericalSuperficalArea_ = nonSphericalModelDict.get<scalar>("nonsphericalParticleSuperficialArea");
        
        //- if solid particle is nonspherical
        scalar Rzero = pow((6.0*nonSphericalVolume_/constant::mathematical::pi),1.0/3.0)/2.0;
        scalar Azero = constant::mathematical::pi*pow((6.0*nonSphericalVolume_/constant::mathematical::pi),2.0/3.0);
        scalar PHI = Azero/nonSphericalSuperficalArea_;
        effD_ = 2.0*Rzero*pow(PHI,-1.0/2.0);
        effRho_ = rho_*pow(PHI,3.0/2.0);
    }
    
    massSphere_ = constant::mathematical::pi*effD_*effD_*effD_*effRho_/6.0;
    
    const dictionary& conductiveHeatTranferDict
    (
        dict.subDict("conductiveHeatTranferProPerties")
    );
    k_ = conductiveHeatTranferDict.get<scalar>("thermalConductivity");
    poissonRatio_ = conductiveHeatTranferDict.get<scalar>("poissonRatio");
    elasticModuli_ = conductiveHeatTranferDict.get<scalar>("elasticModuli");
    
    const dictionary& phaseChangeDict
    (
        dict.subDict("phaseChangeModelProPerties")
    );
    Tm_ = phaseChangeDict.get<scalar>("equilibriumMeltingTemperature");
    Tf_ = phaseChangeDict.get<scalar>("nucleationTemperature");
    
}
	


inline Foam::solidParticleCoupling::solidParticleCoupling
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const vector& U,
    const vector& omega,
    const vector& UCorrect,
    const vector& F,
    const scalar D,
    const scalar T,
    const scalar RWF,
    const scalar CzRatio,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label numSeq,
    const label typeID,
    const label phaseState,
    const label newSolidParticle
)
:
    particle(mesh, coordinates, cellI, tetFaceI, tetPtI),
    U_(U),
    omega_(omega),
    UCorrect_(UCorrect),
    F_(F),
    D_(D),
    T_(T),
    RWF_(RWF),
    CzRatio_(CzRatio),
    numSeq_(numSeq),
    typeID_(typeID),
    phaseState_(phaseState),
    newSolidParticle_(newSolidParticle)
{}

inline Foam::solidParticleCoupling::solidParticleCoupling
(
    const polyMesh& mesh,
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
:
    particle(mesh, position, cellI),
    U_(U),
    omega_(omega),
    UCorrect_(UCorrect),
    F_(F),
    D_(D),
    T_(T),
    RWF_(RWF),
    CzRatio_(CzRatio),
    numSeq_(numSeq),
    typeID_(typeID),
    phaseState_(phaseState),
    newSolidParticle_(newSolidParticle)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//- constantProperties
inline Foam::scalar Foam::solidParticleCoupling::constantProperties::d() const
{
    
    return effD_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::epsilon() const
{
    return epsilon_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::alpha() const
{
    return alpha_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::rho() const
{
    return effRho_;
}


inline Foam::scalar Foam::solidParticleCoupling::constantProperties::Cp() const
{
    return Cp_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::mu() const
{
    return mu_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::tau() const
{
    return tau_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::k() const
{
    return k_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::poissonRatio() const
{
    return poissonRatio_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::elasticModuli() const
{
    return elasticModuli_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::Tm() const
{
    return Tm_;
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::Tf() const
{
    return Tf_;
}

inline bool Foam::solidParticleCoupling::constantProperties::nonSperical() const
{
    return nonSperical_;
}

 
inline Foam::scalar Foam::solidParticleCoupling::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*(0.5*effD_)*(0.5*effD_);
}

inline Foam::scalar Foam::solidParticleCoupling::constantProperties::massSphere() const
{
    return massSphere_;
}


// - particle constructor access
inline const Foam::vector& Foam::solidParticleCoupling::U() const
{
    return U_;
}

inline const Foam::vector& Foam::solidParticleCoupling::omega() const
{
    return omega_;
}

inline const Foam::vector& Foam::solidParticleCoupling::UCorrect() const
{
    return UCorrect_;
}

inline const Foam::vector& Foam::solidParticleCoupling::F() const
{
    return F_;
}

inline Foam::scalar Foam::solidParticleCoupling::D() const
{
    return D_;
}

inline Foam::scalar Foam::solidParticleCoupling::T() const
{
    return T_;
}

inline Foam::scalar Foam::solidParticleCoupling::RWF() const
{
    return RWF_;
}

inline Foam::scalar Foam::solidParticleCoupling::CzRatio() const
{
    if(CzRatio_>1.0 || CzRatio_<0.0)
    {
        FatalErrorIn("solidParticleCoupling::CzRatio()")
        << "    Warning! CzRatio is out of range! "
        << "    It should be in the range between 0 and 1. "
        << "    Please check!"
        << nl
        << exit(FatalError); 
    }
    
    return CzRatio_;
}

inline Foam::label Foam::solidParticleCoupling::numSeq() const
{
    return numSeq_;
}

inline Foam::label Foam::solidParticleCoupling::typeID() const
{
    return typeID_;
}

inline Foam::label Foam::solidParticleCoupling::phaseState() const
{
    return phaseState_;
}

inline Foam::label Foam::solidParticleCoupling::newSolidParticle() const
{
    return newSolidParticle_;
}

//- particle constructor Edit
inline Foam::vector& Foam::solidParticleCoupling::U()
{
    return U_;
}

inline Foam::vector& Foam::solidParticleCoupling::omega()
{
    return omega_;
}

inline Foam::vector& Foam::solidParticleCoupling::UCorrect()
{
    return UCorrect_;
}

inline Foam::vector& Foam::solidParticleCoupling::F()
{
    return F_;
}

inline Foam::scalar& Foam::solidParticleCoupling::D()
{
    return D_;
}

inline Foam::scalar& Foam::solidParticleCoupling::T()
{
    return T_;
}

inline Foam::scalar& Foam::solidParticleCoupling::RWF()
{
    return RWF_;
}

inline Foam::scalar& Foam::solidParticleCoupling::CzRatio()
{
    return CzRatio_;
}

inline Foam::label& Foam::solidParticleCoupling::numSeq()
{
    return numSeq_;
}

inline Foam::label& Foam::solidParticleCoupling::typeID()
{
    return typeID_;
}

inline Foam::label& Foam::solidParticleCoupling::phaseState()
{
    return phaseState_;
}

inline Foam::label& Foam::solidParticleCoupling::newSolidParticle()
{
    return newSolidParticle_;
}

inline Foam::solidParticleCoupling::trackingData::trackPart Foam::solidParticleCoupling::trackingData::part() const
{
    return part_;
}

inline Foam::solidParticleCoupling::trackingData::trackPart& Foam::solidParticleCoupling::trackingData::part()
{
    return part_;
}

// ************************************************************************* //
