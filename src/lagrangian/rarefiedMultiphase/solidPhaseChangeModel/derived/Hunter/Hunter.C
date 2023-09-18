/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "Hunter.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Hunter, 0);

addToRunTimeSelectionTable
(solidPhaseChangeModel, Hunter, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Hunter::Hunter
(
//     const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidPhaseChangeModel(spc, dict)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Hunter::active() const
{
    return true;
}


void Hunter::initialConfiguration()
{

}

void Hunter::temperatureCorrection
(
    solidParticleCoupling& pSolid,
    const scalar& QdeltaC,
    const scalar& deltaT
)
{
    //- update the particle phase
    phaseCheck(pSolid);
    
    label typeIdSolid = pSolid.typeID();
    //- solid particle mass
    scalar massSolid = spc_.constSolidProps(typeIdSolid).massSphere();
    
    scalar deltaTemperature = 0.0;
    
    if(pSolid.phaseState()==1)
    {
        //- process of solidification
        
        scalar deltaR1 = -spc_.mesh().time().deltaTValue()
                        *spc_.phaseChangeModelConst(spc_.materialList()[typeIdSolid])
                        *pow((spc_.constSolidProps(typeIdSolid).Tm()-pSolid.T()),1.8)
                        /(0.5*spc_.constSolidProps(typeIdSolid).d());
        
//         Info<<"deltaR1 = "<<deltaR1<<endl;
                        
                        
        scalar r1New = pSolid.CzRatio() + deltaR1;
        
        scalar deltaCuber1 = pow(r1New,3.0)-pow(pSolid.CzRatio(),3.0);
        
        //- latent heat of fusion hf=1.07e6 J/kg
        deltaTemperature = QdeltaC * deltaT / (spc_.constSolidProps(typeIdSolid).Cp() * massSolid) - (spc_.particleLatentHeatOfFusion(spc_.materialList()[typeIdSolid])/spc_.constSolidProps(typeIdSolid).Cp())*deltaCuber1;
        
        
        scalar Nn = 50.0*deltaTemperature/spc_.constSolidProps(typeIdSolid).Tm();
        
        //- to get the smallest integer for N > 50*deltaTsolid/Tm
        label N(Nn+1);
        
//         Info<<"N ="<<N<<endl;
        if(N > 1)
        {
            //- if N > 1, the initial values of Tsolid and r1 are reassigned to the particle
            //- but here, pSolid.Tsolid() and pSolid.CzRatio have not been changed yet.
            scalar deltaTlocal = spc_.mesh().time().deltaTValue()/N;
        
            for(label i=0; i < N; ++i)
            {
                deltaR1 = -deltaTlocal
                            *spc_.phaseChangeModelConst(spc_.materialList()[typeIdSolid])
                            *pow((spc_.constSolidProps(typeIdSolid).Tm()-pSolid.T()),1.8)
                            /(0.5*spc_.constSolidProps(typeIdSolid).d());
                
//                 Info<<"deltaR1Inside = "<<deltaR1<<endl;
                            
                r1New = pSolid.CzRatio() + deltaR1;
                
                deltaCuber1 = pow(r1New,3.0)-pow(pSolid.CzRatio(),3.0);
                
                deltaTemperature = QdeltaC * deltaTlocal / (spc_.constSolidProps(typeIdSolid).Cp() * massSolid) - (1.07e6/spc_.constSolidProps(typeIdSolid).Cp())*deltaCuber1;
                
//                 Info<<"deltaTemperatureInside = "<<deltaTemperature<<endl;
                //- assignment of particle temperature and CzRatio
                pSolid.T() += deltaTemperature;
                pSolid.CzRatio() = r1New;
            }
        }
        else
        {
            pSolid.T() += deltaTemperature;
            
            pSolid.CzRatio() = r1New;
        }
        
    }
    else if(pSolid.phaseState()==2 && QdeltaC > 0.0)
    {
        //- process of melting
        
        //- update the CzRatio through Equation (12) in "Monte Carlo Simulation of a Rarefied Multiphase Plume Flow"
        //- with deltaT equaling to zero
        pSolid.CzRatio() = pow((QdeltaC*deltaT/(massSolid*(1.07e6))+pow((pSolid.CzRatio()),3.0)),1.0/3.0);
        
        if(pSolid.CzRatio() > 1)
        {
            deltaTemperature = (1.07e6/spc_.constSolidProps(typeIdSolid).Cp())*(pow(pSolid.CzRatio(),3.0)-1.0);
            
            pSolid.CzRatio() = 1.0;
        }
        else
        {
            deltaTemperature = 0.0;
        }
        
        pSolid.T() += deltaTemperature;
        
    }
    else
    {
        //- neither solidifying nor melting
        pSolid.T() += (QdeltaC * deltaT / (spc_.constSolidProps(typeIdSolid).Cp() * massSolid));
    }
    
    //- update the particle phase
    phaseCheck(pSolid);
    
    spc_.particleSizeCorrection(pSolid);
}

void Hunter::phaseCheck(solidParticleCoupling& pSolid)
{
    label typeIdSolid = pSolid.typeID();
    
    if(pSolid.T()< spc_.constSolidProps(typeIdSolid).Tf())
    {
        //- when solid particle temperature is under nucleation temperature
        //- simply considered two cases:
        //- 1. it is supercooled and in solidification process 
        //- 2. It is in pure solid phase
        
        if(pSolid.CzRatio() ==1.0 || pSolid.CzRatio() > 0.0)
        {
            //- the core is liquid, unsteady state
            pSolid.phaseState() = 1;
        }
        else if(pSolid.CzRatio() < 0.0)
        {
            //-reset CzRatio
            pSolid.CzRatio() = 0.0;
            //- pure solid phase
            pSolid.phaseState() = 0;
        }
        else
        {
            //- pure solid phase
            pSolid.phaseState() = 0;
        }
    }
    else if(pSolid.T()< spc_.constSolidProps(typeIdSolid).Tm())
    {
        if(pSolid.CzRatio()<1.0 && pSolid.CzRatio()>0.0)
        {
            //- the core is liquid, unsteady state
            pSolid.phaseState() = 1;
        }
        else if(pSolid.CzRatio() < 0.0)
        {
            //-reset CzRatio
            pSolid.CzRatio() = 0.0;
            //- pure solid phase
            pSolid.phaseState() = 0;
        }
    }
    else if(pSolid.T() >= spc_.constSolidProps(typeIdSolid).Tm())
    {
        //- when solid particle temperature is over melting temperature
        //- simply considered two cases:
        //- 1. it is melting 2. It is in pure liquid phase
        if(pSolid.CzRatio() >= 0.0 && pSolid.CzRatio() < 1.0)
        {
            //- the core is solid, unsteady state
            //- in the case of pure solid heated over melting temperature
            //- means it starts to melt,CzRatio==0
            pSolid.phaseState() = 2;
        }
        else
        {
            //- pure liquid phase
            pSolid.phaseState() = 3;
        }
    }
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
