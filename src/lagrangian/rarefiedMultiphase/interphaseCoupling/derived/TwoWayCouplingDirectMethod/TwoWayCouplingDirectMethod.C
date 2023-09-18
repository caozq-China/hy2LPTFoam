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

#include "TwoWayCouplingDirectMethod.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(TwoWayCouplingDirectMethod, 0);
    addToRunTimeSelectionTable(InterphaseCoupling, TwoWayCouplingDirectMethod, dictionary);



Foam::TwoWayCouplingDirectMethod::TwoWayCouplingDirectMethod
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    InterphaseCoupling(spc, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::TwoWayCouplingDirectMethod::monatomicForceTransferToSolidParticle
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

Foam::scalar Foam::TwoWayCouplingDirectMethod::monatomicEnergyTransferToSolidParticle
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

Foam::vector Foam::TwoWayCouplingDirectMethod::polyatomicForceTransferToSolidParticle
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

Foam::scalar Foam::TwoWayCouplingDirectMethod::polyatomicEnergyTransferToSolidParticle
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

void Foam::TwoWayCouplingDirectMethod::moleculePostCollisionVelocityUpdate
(
    const solidParticleCoupling& pP,
    dsmcParcel& pG
)
{
    //label typeIdG = pG.typeId();
    const label typeIdSolid = pP.typeID();
    
    const label typeIdGas = pG.typeId();
  
    const scalar mDsmc = spc_.dsmcCloudReference()->constProps(typeIdGas).mass();
    
    vector& UM = pG.U();
    const vector& UP = pP.U();
    
    const scalar& Tsolid = pP.T();

    //- Effect on particles

    //vector CR(100,10,20);
    const vector CR = UM - UP;
    
    scalar cR = mag(CR);

    scalar uR = CR.x();
    scalar vR = CR.y();
    scalar wR = CR.z();
    
    scalar cRstar;
    
    vector CRstar;
    
    const scalar tauSolid =  spc_.constSolidProps(typeIdSolid).tau();    

    // - Effect on gas molecules
    // - "Development of a Two-way Coupled Model for Two-phase Rarefied Flow"
    // - Jonathan M. Burt and Iain D. Boyd    5-8 January 2004
    //Random rndGen_typeCollision(clock::getTime());
    // azimuthal angle
    scalar epsilonAngle;

    //- determine deflection angle
    scalar deltaAngle;
    
    if(spc_.rndGenS().sample01<scalar>() <= 1.0 - tauSolid)
    {
        // - SPECULAR REFLECTION
       cRstar = cR;
        
       epsilonAngle = spc_.rndGenS().sample01<scalar>() * 2.0 * pi;
       
       scalar cosElevationAngle = 2.0*spc_.rndGenS().sample01<scalar>()-1;
       scalar sinElevationAngle = sqrt(1-sqr(cosElevationAngle));
       
       CRstar.x() = cosElevationAngle*cRstar;
       CRstar.y() = sinElevationAngle*cos(epsilonAngle)*cRstar;
       CRstar.z() = sinElevationAngle*sin(epsilonAngle)*cRstar;
       
       //Info<<"epsilonAngle = "<< epsilonAngle * 180/ pi << endl;
    }
    else 
    {
        // - DIFFUSE REFLECTION
        //0.72269 is the maximum of SixthOrderPolynomial cConst >=pdf/distribution function of sample
        //scalar cConst = 0.72269*pi;
        
        scalar SixthOrderPolynomial;
        do
        {
            
            deltaAngle = spc_.rndGenS().sample01<scalar>() * pi;//generated sample
            
            SixthOrderPolynomial = 0.02042*pow(deltaAngle,6)-0.2515*pow(deltaAngle,5)+
                                            1.104*pow(deltaAngle,4)-1.903*pow(deltaAngle,3)
                                            +0.4938*pow(deltaAngle,2)+1.248*deltaAngle;//pdf
                                            
        } while (spc_.rndGenS().sample01<scalar>() > SixthOrderPolynomial/0.72269);

        epsilonAngle = spc_.rndGenS().sample01<scalar>() * (2.0 * pi);
       

        scalar cRstarPDF;
        
        const scalar beta = sqrt(mDsmc/(2.0*physicoChemical::k.value()*Tsolid));
        
        const scalar CRstarPDFmax = 3.0*beta*sqrt(1.5)*exp(-1.5);
            
        do
        {
            //- the sample equation of cRstar is not mentioned in Burt and Boyd's paper
            //- this equation is derived from 
            //- "Simulation of rocket plume and lunar dust using DSMC method"
            //cRstar = sqrt(-2.0*log(rndGenS_.scalar01()))/beta;
            
            cRstar = 4.0*spc_.rndGenS().sample01<scalar>()/ beta;
            
            //- equation 12 in Burt and Boyd's paper
            cRstarPDF = 2.0*pow(beta,4.0)*pow(cRstar,3.0)*exp(-sqr(beta*cRstar));
                                                                            
        } while (spc_.rndGenS().sample01<scalar>() > cRstarPDF/CRstarPDFmax);
        
        //Info<<"cRstar = "<< cRstar<<endl;
    
        // determine  gas post-collision relative velocity vector component
        //- correct version (There are typos in the velocity update formula in "Development of 
        //- a Two-way Coupled Model for Two-phase Rarefied Flow")
        scalar uR_postCollisionRelU = (cRstar/cR)*(-uR*cos(deltaAngle)-sin(deltaAngle)
                                            *sin(epsilonAngle)*sqrt(vR*vR + wR*wR));

        scalar vR_postCollisionRelU = (cRstar/cR)*(-vR*cos(deltaAngle)-sin(deltaAngle)*
                                    (cR*wR*cos(epsilonAngle)-uR*vR*sin(epsilonAngle))/sqrt(vR*vR + wR*wR));

        scalar wR_postCollisionRelU = (cRstar/cR)*(-wR*cos(deltaAngle)+sin(deltaAngle)*
                                    (cR*vR*cos(epsilonAngle)+uR*wR*sin(epsilonAngle))/sqrt(vR*vR + wR*wR));
                                    
        /*
        //- version of Burt's PhD thesis
        scalar uR_postCollisionRelU = (cRstar/cR)*(-uR*cos(deltaAngle)+sin(deltaAngle)
                                            *sin(epsilonAngle)*sqrt(vR*vR + wR*wR));

        scalar vR_postCollisionRelU = (cRstar/cR)*(-vR*cos(deltaAngle)-sin(deltaAngle)*
                                    (cR*wR*cos(epsilonAngle)+uR*vR*sin(epsilonAngle))/sqrt(vR*vR + wR*wR));

        scalar wR_postCollisionRelU = (cRstar/cR)*(-wR*cos(deltaAngle)+sin(deltaAngle)*
                                    (cR*vR*cos(epsilonAngle)-uR*wR*sin(epsilonAngle))/sqrt(vR*vR + wR*wR));
	*/

        CRstar.x() = uR_postCollisionRelU;
        CRstar.y() = vR_postCollisionRelU;
        CRstar.z() = wR_postCollisionRelU;
    
        //- for a diffuse reflection, the molecule rotational energy should be updated
        if(spc_.dsmcCloudReference()->constProps(pG.typeId()).rotationalDoF() > 0)
        {
          scalar randNum;
          do
          {
              randNum = spc_.rndGenS().sample01<scalar>();
              
          }while(randNum == 0.0);
          
          pG.ERot() = -log(randNum)*physicoChemical::k.value()*Tsolid;
        }
    }
    
    //Info<<"mag(CRstar) = "<< mag(CRstar)<<endl;
    // - update the gas molecule post collision velocity
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

    ////- resign post collision velocity back to atoms/molecules parcel
    pG.U() = UM;
}


bool Foam::TwoWayCouplingDirectMethod::enableMoleculeVelocityUpdate()
{
    return true;
}

}//end of Foam

// ************************************************************************* //
