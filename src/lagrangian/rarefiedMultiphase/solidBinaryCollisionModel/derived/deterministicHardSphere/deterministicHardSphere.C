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

#include "deterministicHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(deterministicHardSphere, 0);
    addToRunTimeSelectionTable(solidBinaryCollisionModel, deterministicHardSphere, dictionary);
};



Foam::deterministicHardSphere::deterministicHardSphere
(
    const dictionary& dict,
    solidParticleCouplingCloud& spc
)
:
    solidBinaryCollisionModel(dict, spc),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    eN_(coeffDict_.get<scalar>("CoeffNormalResituation")),
    eT_(coeffDict_.get<scalar>("CoeffTangentialResituation")),
    f_(coeffDict_.get<scalar>("CoeffFriction")),
    kn_(coeffDict_.get<scalar>("normalSpringConst")),
    etaN_(coeffDict_.get<scalar>("CoeffDamping"))
    //enableParticleSlide_(coeffDict_.get<Switch>("enableParticleSlide"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::deterministicHardSphere::active() const
{
    return true;
}


void Foam::deterministicHardSphere::collide
(
    solidParticleCoupling& pI,
    solidParticleCoupling& pJ
)
{

    label typeIdi = pI.typeID();
    label typeIdj = pJ.typeID();
    vector& Ui = pI.U();
    vector& Uj = pJ.U();

    vector& omegaI = pI.omega();
    vector& omegaJ = pJ.omega();
    
//     scalar collisionDiameter_pI = 0;
//     scalar collisionDiameter_pJ = 0;
    //- normal vector
    vector n = (pI.position()-pJ.position())/mag(pI.position()-pJ.position());

    vector rOmegaIJ;
    if(spc_.dsmcCloudReference()->axisymmetric())
    {
        rOmegaIJ=(pI.D()*pow(pI.RWF(),1/3)/2)*omegaI+(pJ.D()*pow(pJ.RWF(),1/3)/2)*omegaJ;
//         collisionDiameter_pI = pI.D()*pow((pI.RWF()*spc_.nSolidParticles()),1/3);
//         collisionDiameter_pJ = pJ.D()*pow((pJ.RWF()*spc_.nSolidParticles()),1/3);

        //- restitution coefficient of a parcel
        eN_ = pow(eN_,sqrt(pI.RWF()*spc_.nSolidParticles()));
        eT_ = pow(eT_,sqrt(pI.RWF()*spc_.nSolidParticles()));

    }
    else
    {
        rOmegaIJ=(pI.D()/2)*omegaI+(pJ.D()/2)*omegaJ;
//         collisionDiameter_pI = pI.D()*pow(spc_.nSolidParticles(),1/3);
//         collisionDiameter_pJ = pJ.D()*pow(spc_.nSolidParticles(),1/3);

        //- restitution coefficient of a parcel
        eN_ = pow(eN_,sqrt(spc_.nSolidParticles()));
        eT_ = pow(eT_,sqrt(spc_.nSolidParticles()));
    }

    //- relative velocity


    vector Vc = Uj - Ui- (rOmegaIJ^n);
    vector Vn = mag(Vc)*n;
    vector Vt = Vc - Vn;

    scalar mI = spc_.constSolidProps(typeIdi).massSphere();
    scalar mJ = spc_.constSolidProps(typeIdj).massSphere();
    
    const scalar mE = mI*mJ/(mI+mJ);




    //- variation of momentum
    vector deltaP = -mE*(1+eN_)*Vn-(2*mE*(1+eT_)*Vt)/7;


    scalar Jt_a = 2*mE*(1+eT_)*mag(Vt)/7;
    scalar Jt_b = f_*mE*(1+eN_)*mag(Vn);
    vector Jt;
    if(Jt_a>Jt_b)
    {
        Jt = -(2/7)*mE*(1+eT_)*Vt;
    }
    else
    {
        Jt = -f_*mE*(1+eN_)*Vn;
    }

    Ui -= deltaP/mI;
    Uj += deltaP/mJ;

    //- "^" means vector cross product
    omegaI -= pI.D()*pow(pI.RWF(),1/3)*(n^Jt)/(2*mI*pI.D()*pow(pI.RWF(),1/3)*pI.D()*pow(pI.RWF(),1/3)/10);
    omegaJ -= pJ.D()*pow(pJ.RWF(),1/3)*(n^Jt)/(2*mJ*pJ.D()*pow(pJ.RWF(),1/3)*pJ.D()*pow(pJ.RWF(),1/3)/10);

//     scalar newKineticEP = 0.5*(pP.U()&pP.U());
//     scalar newKineticEQ = 0.5*(pQ.U()&pQ.U());
//
//     if(e_<1)
//     {
//         pP.T() += (newKineticEP-initialKineticEP)/spc_.constSolidProps(typeIdP).Cp();
//         pQ.T() += (newKineticEQ-initialKineticEQ)/spc_.constSolidProps(typeIdQ).Cp();
//     }
    
}

void Foam::deterministicHardSphere::velocityCorrection
(
    int correctionStep,
    solidParticleCoupling& pI,
    solidParticleCoupling& pJ
)
{
    vector rij=pI.position()-pJ.position();

    vector n = (pI.position()-pJ.position())/mag(pI.position()-pJ.position());

//     scalar collisionDiameter_pI = 0;
//     scalar collisionDiameter_pJ = 0;
    vector rOmegaIJ;
    if(spc_.dsmcCloudReference()->axisymmetric())
    {
        rOmegaIJ=(pI.D()*pow(pI.RWF(),1/3)/2)*pI.omega()+(pJ.D()*pow(pJ.RWF(),1/3)/2)*pJ.omega();
//         collisionDiameter_pI = pI.D()*pow((pI.RWF()*spc_.nSolidParticles()),1/3);
//         collisionDiameter_pJ = pJ.D()*pow((pJ.RWF()*spc_.nSolidParticles()),1/3);
    }
    else
    {
        rOmegaIJ=(pI.D()/2)*pI.omega()+(pJ.D()/2)*pJ.omega();
//         collisionDiameter_pI = pI.D()*pow(spc_.nSolidParticles(),1/3);
//         collisionDiameter_pJ = pJ.D()*pow(spc_.nSolidParticles(),1/3);
    }

    vector Vc = pJ.U() - pI.U()- (rOmegaIJ^n);

    label typeIdi = pI.typeID();
//     label typeIdj = pJ.typeID();

    scalar mI = spc_.constSolidProps(typeIdi).massSphere();
//     scalar mJ = spc_.constSolidProps(typeIdj).massSphere();

    if(correctionStep==0)
    {
        Vc=Zero;
    }

    vector Vn = mag(Vc)*n;

    scalar deltan = (pI.D()*pow(pI.RWF(),1/3)+pJ.D()*pow(pJ.RWF(),1/3))/2 - mag(rij);//- overlap length
    if((rij&Vn)>0 && (deltan>0))
    {
        scalar magVc = mag(Vc);
        if(correctionStep==0)
        {
            magVc = -mag(Vn);
        }

        //- time step
        const scalar deltaT = spc_.dsmcCloudReference()->mesh().time().deltaTValue();

        scalar koverlap=deltaT*kn_/mI;
        scalar kvelocity=deltaT*etaN_/mI;

        scalar deltaV = -koverlap*deltan+kvelocity*(magVc-mag(Vn));
        pI.U() -= n*deltaV;
        pJ.U() += n*deltaV;
//         noMoreCollide = false;
    }
}

const Foam::dictionary& Foam::deterministicHardSphere::coeffDict() const
{
    return coeffDict_;
}
// ************************************************************************* //
