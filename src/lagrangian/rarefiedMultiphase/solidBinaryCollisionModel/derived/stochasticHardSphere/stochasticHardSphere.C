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

#include "stochasticHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(stochasticHardSphere, 0);
    addToRunTimeSelectionTable(solidBinaryCollisionModel, stochasticHardSphere, dictionary);
};



Foam::stochasticHardSphere::stochasticHardSphere
(
    const dictionary& dict,
    solidParticleCouplingCloud& spc
)
:
    solidBinaryCollisionModel(dict, spc),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    e_(coeffDict_.get<scalar>("CoeffResituation"))
    //f_(coeffDict_.get<scalar>("CoeffFriction")),
    //enableParticleSlide_(coeffDict_.get<Switch>("enableParticleSlide"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::stochasticHardSphere::active() const
{
    return true;
}


void Foam::stochasticHardSphere::collide
(
    solidParticleCoupling& pP,
    solidParticleCoupling& pQ
)
{

    label typeIdP = pP.typeID();
    label typeIdQ = pQ.typeID();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
    scalar mP = spc_.constSolidProps(typeIdP).massSphere();
    scalar mQ = spc_.constSolidProps(typeIdQ).massSphere();
    
    const vector cM = (mP*UP+mQ*UQ)/(mP+mQ);
    
    //- generate a unit vector on a sphere
    scalar cosTheta = 2.0*spc_.rndGenS().sample01<scalar>() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*spc_.rndGenS().sample01<scalar>();
    
    vector vecE(cosTheta, sinTheta*cos(phi),sinTheta*sin(phi));
    
    pP.U() = cM + e_*mag(UP - cM)*vecE;
    pQ.U() = cM - e_*mag(UQ - cM)*vecE;
    
    scalar initialKineticEP = 0.5*(UP&UP);
    scalar initialKineticEQ = 0.5*(UQ&UQ);
    
 /* 
 
 //- The code below is based on Chapter 5 of "Multiphase Flows with Droplets and Particles (2nd Edition) Crowe C, Sommerfeld M, et.al".
 //- And they are based on deterministic collision detection method, which is not suitable for stochastic collision detection algorithm, e.g. no time counter.
 
    vector G0 = UP - UQ;
    
    //- unit normal vector
    vector n = (pP.position() - pQ.position())/mag((pP.position() - pQ.position()));
    
    vector G0ct = G0 - (G0 & n) * n;
    
    //- unit tangential vector
    vector t = G0ct/mag(G0ct);
    
    scalar mP = spc_.constSolidProps(typeIdP).massSphere();
    scalar mQ = spc_.constSolidProps(typeIdQ).massSphere();
    
//     vector cM = (mP*UP+mQ*UQ)/(mP+mQ);
//     
//     //-generate a unit vector uniformly distributed on a sphere
//     scalar theta = 2.0*pi*spc_.rndGenS().scalar01();
//     //- generate a number between [-1,1]
//     scalar Z = spc_.rndGenS().scalar01()*2.0-1.0;
//     vector xn(sqrt(1-sqr(Z))*cos(theta),sqrt(1-sqr(Z))*sin(theta),Z);
    
    
    if(enableParticleSlide_)
    {
         pP.U() = UP - (n + f_*t)*(n & G0)*(1.0+e_)*mQ/(mP+mQ);
         pQ.U() = UQ + (n + f_*t)*(n & G0)*(1.0+e_)*mP/(mP+mQ);
    }
    else
    {
         pP.U() = UP - ((1.0+e_)*(n & G0)*n+2.0*mag(G0ct)*t/7.0)*mQ/(mP+mQ);
         pQ.U() = UQ + ((1.0+e_)*(n & G0)*n+2.0*mag(G0ct)*t/7.0)*mP/(mP+mQ);
//            pP.U() = cM + e_*mag(UP-cM)*xn;
//            pQ.U() = cM - e_*mag(UP-cM)*xn;
            
    }
*/    
    scalar newKineticEP = 0.5*(pP.U()&pP.U());
    scalar newKineticEQ = 0.5*(pQ.U()&pQ.U());
    
    if(e_<1)
    {
        pP.T() += (newKineticEP-initialKineticEP)/spc_.constSolidProps(typeIdP).Cp();
        pQ.T() += (newKineticEQ-initialKineticEQ)/spc_.constSolidProps(typeIdQ).Cp();
    }
    
}

void Foam::stochasticHardSphere::velocityCorrection
(
    int nCorrectionStep,
    solidParticleCoupling& pP,
    solidParticleCoupling& pQ
)
{
}


const Foam::dictionary& Foam::stochasticHardSphere::coeffDict() const
{
    return coeffDict_;
}
// ************************************************************************* //
