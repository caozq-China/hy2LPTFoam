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

#include "TwoWayCouplingIndirectMethod.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(TwoWayCouplingIndirectMethod, 0);
    addToRunTimeSelectionTable(InterphaseCoupling, TwoWayCouplingIndirectMethod, dictionary);



Foam::TwoWayCouplingIndirectMethod::TwoWayCouplingIndirectMethod
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    InterphaseCoupling(spc, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::TwoWayCouplingIndirectMethod::monatomicForceTransferToSolidParticle
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

Foam::scalar Foam::TwoWayCouplingIndirectMethod::monatomicEnergyTransferToSolidParticle
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

Foam::vector Foam::TwoWayCouplingIndirectMethod::polyatomicForceTransferToSolidParticle
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

Foam::scalar Foam::TwoWayCouplingIndirectMethod::polyatomicEnergyTransferToSolidParticle
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

void Foam::TwoWayCouplingIndirectMethod::moleculePostCollisionVelocityUpdate
(
    const solidParticleCoupling& pP,
    dsmcParcel& pG
)
{
//     const label typeIdP = pP.typeIDsolid();
 
//     label typeIdGas = pG.typeId();

    const scalar mDsmc = spc_.dsmcCloudReference()->constProps(pG.typeId()).mass();

    vector& UM = pG.U();

    const vector& UP = pP.U();
    
    const scalar dsolid = pP.D();
    
    vector CR = UM - UP;
    //vector CR(100,10,20);
    vector CRstar;
    
    vector CRunit = CR / mag(CR);
    const vector XminusUnit(-1.0,0.0,0.0);
    const vector Xunit(1.0,0.0,0.0);
    //rotation matrix T1= rotationTensor(VectorbeforeTransformation,vectorAfterTransformation)
    tensor T1 = rotationTensor(CRunit,XminusUnit);
    tensor T1Inv = inv(T1);
    
    //impact parameter
    scalar b = dsolid * sqrt(spc_.rndGenS().sample01<scalar>()) / 2.0;
    
    scalar phi = 2.0 * pi * spc_.rndGenS().sample01<scalar>();
    
    //reflection point
    vector OM(sqrt(sqr(dsolid / 2.0) - sqr(b)),b * cos(phi),b * sin(phi));
    vector OMunit = OM / mag(OM);
    
    tensor T2 = rotationTensor(OMunit,Xunit);
    tensor T2inv = inv(T2);
    
    //- p is unit vector normalized from vector p
    //in "Simulation of rocket plume and lunar dust using DSMC method"
    // by "He Xiaoying, He Bijiao, Cai Guobiao" Acta Astronautica 70(2012) 100-111
    vector p = T2 & XminusUnit;
    
    // pPost is the post collision velocity in the normal coordinate syetem
    vector pPost;
    
    scalar tauSolid =  spc_.constSolidProps(pP.typeID()).tau();  
    
    // determine the deflection type
    if(spc_.rndGenS().sample01<scalar>() <= 1.0 - tauSolid)
    {
        //- SPECULAR REFLECTION
        pPost.x() = -p.x() * mag(CR);
        pPost.y() = p.y() * mag(CR);
        pPost.z() = p.z() * mag(CR);
    }
    else 
    {
        //- DIFFUSE REFLECTION
        // - The acceptance-rejection method -//
        //sample of velocity component U
        scalar beta = sqrt(mDsmc/(2.0 *physicoChemical::k.value() * pP.T()));
        //- normal component in normal corrdinate
        scalar Usample;
        //tangential component in normal corrdinate
        scalar Vsample;
        
        //- randomly sample and avoid the mathematical error of log(0)
        scalar randNum1;
        do
        {
            randNum1 = spc_.rndGenS().sample01<scalar>();
            
        }while(randNum1 == 0.0);
        
        Usample = sqrt(-log(randNum1)) / beta;
        
        //- randomly sample and avoid the mathematical error of log(0)
        scalar randNum2;
        do
        {
            randNum2 = spc_.rndGenS().sample01<scalar>();
            
        }while(randNum2 == 0.0);
        
        Vsample = sqrt(-log(randNum2)) / beta;

        //diffusive reflection
        scalar alpha = spc_.rndGenS().sample01<scalar>() * 2.0 * pi;
        pPost.x() = Usample;
        pPost.y() = Vsample*cos(alpha);
        pPost.z() = Vsample*sin(alpha);
        
        //- for a diffuse reflection, the molecule rotational energy should be updated
        if(spc_.dsmcCloudReference()->constProps(pG.typeId()).rotationalDoF() > 0)
        {
            scalar randNum;
            do
            {
                randNum = spc_.rndGenS().sample01<scalar>();
                
            }while(randNum < VSMALL);
            
            pG.ERot() = -log(randNum)*physicoChemical::k.value()*pP.T();
        }
    }
    //- g = T1*Cr
    //- p = T2*g so p = T2*T1*Cr
    // ----> so Cr = T1Inv*T2Inv*p
    
    CRstar = T1Inv &  T2inv & pPost;
    
    
    
    //Info<<"cRstar = "<< mag(CRstar)<<endl;
    //- assign new gas molecule velocity back to the parcel
    UM = CRstar + UP;
    /*
    //- For validation output
    scalar COSpolarAngle = CRstar.y()/mag(CRstar);
    
    //scalar SinPolarAngle = sqrt(1-sqr(COSpolarAngle));
    
    scalar azimuthalAngle;
    
    azimuthalAngle = (atan2(CRstar.z(),CRstar.x())*180/pi);
    if(azimuthalAngle<0)
    {
        azimuthalAngle = azimuthalAngle +360;
    }
    
    
    Info<<"cosPhi = "<< COSpolarAngle<<endl;
    Info<<"azimuthalAngle = "<< azimuthalAngle << endl;
    */
    
    
    //- resign post collision velocity back to atoms/molecules parcel
    pG.U() = UM;
        
        
}

bool Foam::TwoWayCouplingIndirectMethod::enableMoleculeVelocityUpdate()
{
    return true;
}


}//end of Foam

// ************************************************************************* //
