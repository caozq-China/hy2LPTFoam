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
    FITNESS FOR A PARTICULAR PcRPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "OneWayCoupling.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(OneWayCoupling, 0);
    addToRunTimeSelectionTable(InterphaseCoupling, OneWayCoupling, dictionary);



Foam::OneWayCoupling::OneWayCoupling
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    InterphaseCoupling(spc, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::vector Foam::OneWayCoupling::monatomicForceTransferToSolidParticle
(
    const solidParticleCoupling& pP,
    const dsmcParcel& pG,
    const scalar& cellVolume,
    const scalar& gasRWF
) const
{
    //- solid particle properties
    label typeIdSolid = pP.typeID();
    label typeIdGas = pG.typeId();
    scalar mDsmc = spc_.dsmcCloudReference()->constProps(typeIdGas).mass();
    scalar nParticle = spc_.dsmcCloudReference()->nParticle()*gasRWF;
    scalar dsmcParticleMassDensity = mDsmc * nParticle / cellVolume;
    
    scalar epsilonSolid =  spc_.constSolidProps(typeIdSolid).epsilon();
    scalar alphaSolid =  spc_.constSolidProps(typeIdSolid).alpha();
    vector Usolid = pP.U();
    scalar Tsolid =  pP.T();
//     scalar dsolid = spc_.constSolidProps(typeIdSolid).dSolid();
    scalar dsolid = pP.D();
    
    //-relative velocity
    vector Ur = (pG.U() - Usolid);
    //scalar magUr = mag(Ur);
    
    //- force addition
    return dsmcParticleMassDensity*(pi*sqr(dsolid/2))*Ur*
                ((1.0+(4.0/9.0)*(1.0-epsilonSolid)*(1.0-alphaSolid))*
                mag(Ur)+(1.0-epsilonSolid)*alphaSolid*(sqrt(pi)/3.0)*sqrt((2.0*physicoChemical::k.value()*Tsolid)/mDsmc));
}

Foam::scalar Foam::OneWayCoupling::monatomicEnergyTransferToSolidParticle
(
    const solidParticleCoupling& pP,
    const dsmcParcel& pG,
    const scalar& cellVolume,
    const scalar& gasRWF
) const
{
    //- solid particle properties
    label typeIdSolid = pP.typeID();
    label typeIdGas = pG.typeId();
    scalar mDsmc = spc_.dsmcCloudReference()->constProps(typeIdGas).mass();
    scalar nParticle = spc_.dsmcCloudReference()->nParticle()*gasRWF;
    scalar dsmcParticleMassDensity = mDsmc * nParticle / cellVolume;
    
    scalar epsilonSolid =  spc_.constSolidProps(typeIdSolid).epsilon();
    scalar alphaSolid =  spc_.constSolidProps(typeIdSolid).alpha();
    vector Usolid = pP.U();
    scalar Tsolid =  pP.T();
//     scalar dsolid = spc_.constSolidProps(typeIdSolid).dSolid();
    scalar dsolid = pP.D();
    
    //-relative velocity
    vector Ur = (pG.U() - Usolid);
    //scalar magUr = mag(Ur);

    //- energy addition    
    return (1.0-epsilonSolid)*alphaSolid*dsmcParticleMassDensity*
                (pi*sqr(dsolid/2.0))*mag(Ur)*
                (0.5*sqr(mag(Ur))-(2.0*physicoChemical::k.value()*Tsolid)/mDsmc);
}

Foam::vector Foam::OneWayCoupling::polyatomicForceTransferToSolidParticle
(
    const solidParticleCoupling& pP,
    const dsmcParcel& pG,
    const scalar& cellVolume,
    const scalar& gasRWF
) const
{
    //- solid particle properties
    label typeIdSolid = pP.typeID();
    label typeIdGas = pG.typeId();
    scalar mDsmc = spc_.dsmcCloudReference()->constProps(typeIdGas).mass();
    scalar nParticle = spc_.dsmcCloudReference()->nParticle()*gasRWF;
    
    scalar tauSolid =  spc_.constSolidProps(typeIdSolid).tau();    
    vector Usolid = pP.U();
    scalar Tsolid =  pP.T();
//     scalar dsolid = spc_.constSolidProps(typeIdSolid).dSolid();
    scalar dsolid = pP.D();

    //-relative velocity
    vector Ur = (pG.U() - Usolid);
    
    //-force addition
    return ((pi*sqr(0.5 * dsolid)*nParticle)/cellVolume)*(mDsmc*mag(Ur)+tauSolid*sqrt(2.0*pi*mDsmc*physicoChemical::k.value()*Tsolid)/3.0)*Ur;

}

Foam::scalar Foam::OneWayCoupling::polyatomicEnergyTransferToSolidParticle
(
    const solidParticleCoupling& pP,
    const dsmcParcel& pG,
    const scalar& cellVolume,
    const scalar& gasRWF
) const
{
    //- solid particle properties
    label typeIdSolid = pP.typeID();
    label typeIdGas = pG.typeId();
    scalar mDsmc = spc_.dsmcCloudReference()->constProps(typeIdGas).mass();
    scalar nParticle = spc_.dsmcCloudReference()->nParticle()*gasRWF;
    scalar rotationalDegreeOfFreedom = spc_.dsmcCloudReference()->constProps(typeIdGas).rotationalDoF();
    
    scalar tauSolid =  spc_.constSolidProps(typeIdSolid).tau();    
    vector Usolid = pP.U();
    scalar Tsolid =  pP.T();
//     scalar dsolid = spc_.constSolidProps(typeIdSolid).dSolid();
    scalar dsolid = pP.D();

    //-relative velocity
    vector Ur = (pG.U() - Usolid);
    
    //- energy addition
    return (pi*sqr(0.5*dsolid)*tauSolid*nParticle*mag(Ur)/cellVolume)*(0.5*mDsmc*sqr(mag(Ur))+pG.ERot()-(2.0+0.5*rotationalDegreeOfFreedom)*physicoChemical::k.value()*Tsolid);
}

void Foam::OneWayCoupling::moleculePostCollisionVelocityUpdate
(
    const solidParticleCoupling& pP,
    dsmcParcel& pG
)
{
    
}

bool Foam::OneWayCoupling::enableMoleculeVelocityUpdate()
{
    return false;
}


}//end of Foam

// ************************************************************************* //
