/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Thermophoretic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Thermophoretic, 0);
    addToRunTimeSelectionTable
    (
        Force,
        Thermophoretic,
        dictionary
    );
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Thermophoretic::Thermophoretic
(
    solidParcelCloud& cloud,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Force(cloud, mesh, dict, typeName, false),
    TName_(dict.lookupOrDefault<word>("Tt", "Tt")),
    gradTInterpPtr_(nullptr)
{}


Foam::Thermophoretic::Thermophoretic(const Thermophoretic& gf)
:
    Force(gf),
    gradTInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
Foam::Thermophoretic::~Thermophoretic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::Thermophoretic::cacheFields(const bool store)
{
    static word fName("gradT");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volScalarField& Tc = this->mesh().template
                lookupObject<volScalarField>(TName_);

            volVectorField* gradTPtr = new volVectorField
            (
                fName,
                (fvc::grad(Tc))
            );

            gradTPtr->store();
        }

        const volVectorField& gradT = this->mesh().template
            lookupObject<volVectorField>(fName);

        gradTInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                cloud_.solution().interpolationSchemes(),
                gradT
            ).ptr()
        );
    }
    else
    {
        gradTInterpPtr_.clear();

        if (fieldExists)
        {
            const volVectorField& gradT = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(gradT).checkOut();
        }
    }
}

Foam::forceSuSp Foam::Thermophoretic::calcNonCoupled
(
    const solidParcel& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,//- based on relative velocity
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    const scalar magUc = mag(p.Uc());
    if(magUc>ROOTVSMALL)
    {
        const scalar Ct=2.18;
        const scalar cTheta = 1.22;
        const scalar cThetaRatio = (2-cTheta)/cTheta;
        // Info<<"cThetaRatio="<<cThetaRatio<<endl;
        const scalar MachP = mag(p.Uc()-p.U())*p.Mac()/magUc;
        // Info<<"MachP="<<MachP<<endl;
        const scalar KnP = MachP*sqrt(p.gammac()*pi/2)/Re;
        // Info<<"KnP="<<KnP<<endl;
        // Info<<"p.K()="<<p.K()<<endl;
        const scalar kappaStar = p.kappac()/p.K();
        // Info<<"kappaStar="<<kappaStar<<endl;
        scalar A = 0.0;
        if(KnP<=0.01)
        {
            A = 6*pi*sqr(muc)/p.rhoc()*p.d()*cThetaRatio*(kappaStar+2*KnP*Ct)
                        /(1+6*KnP*cThetaRatio)*(1+2*kappaStar+4*cThetaRatio*Ct);
        }
        else
        {
            A = 0.5*pi*sqr(muc)/p.rhoc()*p.d()/(1.15+KnP);
        }
        vector gradT = gradTInterp().interpolate(p.position(), p.currentTetIndices());
        // Info<<"gradT="<<gradT<<endl;
        // Info<<"p.Tc()="<<p.Tc()<<endl;
        value.Su() = -A*gradT/p.Tc();
        // Info<<"value="<<value<<endl;
    }
    return value;
}


// ************************************************************************* //
