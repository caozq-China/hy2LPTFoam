/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "SinghDrag.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SinghDrag, 0);
    addToRunTimeSelectionTable
    (
        Force,
        SinghDrag,
        dictionary
    );
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Reynolds number based on gas velocity
Foam::scalar Foam::SinghDrag::ReInfinity
(
    const scalar d, 
    const vector Uc_, 
    const scalar rhoc, 
    const scalar muc
) const
{
    return rhoc*mag(Uc_)*d/(muc + ROOTVSMALL);
}


// Foam::scalar Foam::SinghDrag::volumeFractionCorrectionFactor(const scalar thetaCCell) const
// {
//     scalar thetaPCell = 1.0 - thetaCCell;
//     return (1.0+2.0*thetaPCell)/sqr(1.0-thetaPCell);
// }


Foam::scalar Foam::SinghDrag::Cd
(
    const scalar Ma,
    const scalar Re,
    const scalar gamma,//- specific heat ratio 
    const scalar Tc,
    const scalar Tp,
    const scalar omega //- viscosity index
) const
{
    scalar Cd = 0;
    if(Ma>ROOTVSMALL)
    {
        scalar Kn=Ma*sqrt(gamma*pi/2)/Re;
        scalar C0 = 0.271616116;//24.0/sqr(9.4);

        scalar Ccd = 0;


        scalar UsOverU=0;
        scalar TsOverT=0;
        scalar MaS=0;
        scalar WTr=0;
        scalar alpha=0;
        scalar C1 = 0;

        if(Ma<=1)
        {
            C1=1;
            alpha=1;
            MaS=Ma;
            UsOverU=1;
            TsOverT=1;
            WTr = pow(Ma,(2*omega))/Re;
        }
        else
        {
            scalar CdM=0.9;
            C1 = (CdM-C0*pow(1+sqr(gamma-1)/(4*gamma),gamma/(gamma-1)))/(1-(1/(0.356*Ma))*((gamma-1)/(gamma+1)));

            alpha = 1/(0.356*(Ma-1)+1);
            MaS = sqrt(((gamma-1)*sqr(Ma)+2)/(2*gamma*sqr(Ma)-(gamma-1)));
            UsOverU=(2+(gamma-1)*sqr(Ma))/((gamma+1)*sqr(Ma));
            TsOverT=((gamma-1)*sqr(Ma)+2)*(2*gamma*sqr(Ma)-(gamma-1))/sqr(((gamma+1)*Ma));
            scalar TpOverTs=Tp/(TsOverT*Tc);
            WTr = (pow(Ma,(2*omega))/Re)*pow((1+TpOverTs),omega);
        }
        scalar thetaS = pow((1+(gamma-1)*sqr(MaS)/2),(gamma/(gamma-1)));
        scalar ReS_tilde = Re*pow(1/(sqr(alpha)*TsOverT),omega)
                            *pow(thetaS,(gamma+1)/(2*gamma)-omega*(gamma-1)/gamma);

        Ccd = C1*(1-alpha*UsOverU)+C0*thetaS*sqr(1+9.4/sqrt(ReS_tilde));
        
        scalar fKnWr = 1/((1+Kn*(2.514+0.8*exp(-0.55/Kn)))*(1+1.27*WTr));
        scalar Br = WTr*(pow(Ma,(2*omega-1))+1)/pow(Ma,(2*omega-1));
        scalar s = Ma*sqrt(gamma/2);
        scalar Cdfm = (1+2*sqr(s))*exp(-sqr(s))/(pow(s,3)*sqrt(pi))
                        + (4*pow(s,4)+4*sqr(s)-1)*Foam::erf(s)/(2*pow(s,4))
                        + 2*sqrt(pi*Tp/Tc)/(3*s);

        Cd = Ccd*fKnWr*(1/(1+pow(Br,1.8)))+Cdfm*pow(Br,1.8)/(1+pow(Br,1.8));
    }

    return Cd;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SinghDrag::SinghDrag
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, true)
    // thetaC_
    // (
    //     this->mesh().template lookupObject<volScalarField>
    //     (
    //         dict.lookup("gasVolumeFractionRef")
    //     )
    // ),
    // volumeFractionCorrection_(Switch(dict.lookup("volumeFractionCorrection")))
{}


Foam::SinghDrag::SinghDrag
(
    const SinghDrag& df
)
:
    Force(df)
    // thetaC_(df.thetaC_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SinghDrag::~SinghDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::forceSuSp Foam::SinghDrag::calcCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
   
    forceSuSp value(Zero, 0.0);
    const scalar Reinfi = ReInfinity(p.d(),p.Uc(),p.rhoc(),muc);
   
    scalar cd = Cd(p.Mac(),Reinfi,p.gammac(),p.Tc(),p.T(),p.omegac());
    // if(volumeFractionCorrection_)
    // {
    //      scalar thetaCCell(thetaC_[p.cell()]);
         // cd *= volumeFractionCorrectionFactor(thetaCCell);
    // }

    value.Sp() = (mass/p.rho())*0.75*cd*Re*muc/sqr(p.d());
    return value;
}

// ************************************************************************* //
