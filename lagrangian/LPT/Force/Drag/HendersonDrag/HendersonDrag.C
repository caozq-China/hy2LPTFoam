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

#include "HendersonDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HendersonDrag, 0);
    addToRunTimeSelectionTable
    (
        Force,
        HendersonDrag,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::HendersonDrag::molecularSpeedRatio(const scalar Map, const scalar gamma) const
{
    return Map*sqrt(0.5*gamma);
}

Foam::scalar Foam::HendersonDrag::TpOverTc(const scalar Tp, const scalar Tc) const
{
    return Tp/Tc;
}

Foam::scalar Foam::HendersonDrag::Cd
(
    const scalar Ma,
    const scalar Re,
    const scalar gamma,//- specific heat ratio 
    const scalar Tc,
    const scalar Tp
) const
{
    scalar Cd=0;
    const scalar S=molecularSpeedRatio(Ma,gamma);
    if(Ma<1)
    {
         Cd = 24/(Re+S*(4.33+(3.65-1.53*TpOverTc(Tp,Tc))/(1+0.353*TpOverTc(Tp,Tc)))*exp(-0.247*Re/S))
              + exp(-0.5*Ma/sqrt(Re))*((4.5+0.38*(0.03*Re+0.48*sqrt(Re)))/(1+0.03*Re+0.48*sqrt(Re))+0.1*sqr(Ma)+0.2*pow(Ma,8))
              + 0.6*S*(1-exp(-Ma/Re));
    }
    else if(Ma>1.75)
    {
        Cd = (0.9+0.34/sqr(Ma)+1.86*sqrt(Ma/Re)*(2+2/sqr(S)+1.058*sqrt(TpOverTc(Tp,Tc))/S-1/pow(S,4)))/(1+1.86*sqrt(Ma/Re));
    }
    else
    {
        const scalar S1=molecularSpeedRatio(1,gamma);
        
        scalar CdMa1=24/(Re+S1*(4.33+(3.65-1.53*TpOverTc(Tp,Tc))/(1+0.353*TpOverTc(Tp,Tc)))*exp(-0.247*Re/S1))
              + exp(-0.5*1/sqrt(Re))*((4.5+0.38*(0.03*Re+0.48*sqrt(Re)))/(1+0.03*Re+0.48*sqrt(Re))+0.1+0.2)
              + 0.6*S1*(1-exp(-1/Re));
        const scalar S1point75=molecularSpeedRatio(1.75,gamma);
        
        scalar CdMa1point75=(0.9+0.34/3.0625+1.86*sqrt(1.75/Re)*(2+2/sqr(S1point75)+1.058*sqrt(TpOverTc(Tp,Tc))/S1point75-1/pow(S,4)))/(1+1.86*sqrt(1.75/Re));

        Cd = CdMa1+4*(Ma-1)*(CdMa1point75-CdMa1)/3;
    }

    return Cd;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::HendersonDrag::HendersonDrag
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, true)
{}


Foam::HendersonDrag::HendersonDrag
(
    const HendersonDrag& df
)
:
    Force(df)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HendersonDrag::~HendersonDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::forceSuSp Foam::HendersonDrag::calcCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    //- Re is based on gas-particle relative velocity
    scalar cd = Cd(p.Mac(),Re,p.gammac(),p.Tc(),p.T());

    value.Sp() = (mass/p.rho())*0.75*cd*Re*muc/sqr(p.d());
    
    return value;
}


// ************************************************************************* //
