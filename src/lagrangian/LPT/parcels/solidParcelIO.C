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

#include "solidParcel.H"
#include "solidParcelCloud.H"
#include "IOstreams.H"
#include "IOField.H"
#include "vectorFieldIOField.H"
#include "labelFieldIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParcel::solidParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    U_(Zero),
    UCorrect_(Zero),
    RWF_(0.0),
    d_(0.0),
    rho_(0.0),
    T_(0.0),
    K_(0.0),
    Cp_(0.0),
    typeId_(0),
    newParcel_(false),
    active_(true),
    Uc_(Zero),
    rhoc_(0.0),
    muc_(0.0),
    Tc_(0.0),
    kappac_(0.0),
    Cpc_(0.0),
    Mac_(0.0),
    gammac_(0.0),
    omegac_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            is >> UCorrect_;
            RWF_ = readScalar(is);
            d_ = readScalar(is);
            rho_ = readScalar(is);
            T_ = readScalar(is);
            K_ = readScalar(is);
            Cp_ = readScalar(is);
            typeId_= readLabel(is);
            newParcel_ = readLabel(is);
            active_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
                +sizeof(UCorrect_)
                +sizeof(RWF_)
                +sizeof(d_)
                +sizeof(rho_)
                +sizeof(T_)
                +sizeof(K_)
                +sizeof(Cp_)
                +sizeof(typeId_)
                +sizeof(newParcel_)
                +sizeof(active_)
            );
        }
    }

    // Check state of Istream
    is.check(FUNCTION_NAME);
}


void Foam::solidParcel::readFields(Cloud<solidParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<vector> UCorrect(c.fieldIOobject("UCorrect", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UCorrect);

    IOField<scalar> RWF(c.fieldIOobject("RWF", IOobject::MUST_READ));
    c.checkFieldIOobject(c, RWF);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rho);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ));
    c.checkFieldIOobject(c, T);

    IOField<scalar> K(c.fieldIOobject("K", IOobject::MUST_READ));
    c.checkFieldIOobject(c, K);

    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Cp);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newParcel);

    IOField<label> active(c.fieldIOobject("active", IOobject::MUST_READ));
    c.checkFieldIOobject(c, active);

    //--------

    label i = 0;
    forAllIter(Cloud<solidParcel>, c, iter)
    {
        solidParcel& p = iter();

        p.U_ = U[i];
        p.UCorrect_ = UCorrect[i];
        p.RWF_ = RWF[i]; 
        p.d_ = d[i]; 
        p.rho_ = rho[i];
        p.T_ = T[i];
        p.K_ = K[i];
        p.Cp_ = Cp[i];
        p.typeId_ = typeId[i];
        
        p.newParcel_ = newParcel[i];
        p.active_ = active[i];

        i++;
    }
}


void Foam::solidParcel::writeFields(const Cloud<solidParcel>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> UCorrect(c.fieldIOobject("UCorrect", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("RWF", IOobject::NO_READ), np);
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> K(c.fieldIOobject("K", IOobject::NO_READ), np);
    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    
    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::NO_READ), np);
    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    
    //-------

    label i = 0;
    forAllConstIter(solidParcelCloud, c, iter)
    {
        const solidParcel& p = iter();
        
        U[i] = p.U();
        UCorrect[i] = p.UCorrect();
        RWF[i] = p.RWF();
        d[i] = p.d();
        rho[i] = p.rho();
        T[i] = p.T();
        K[i] = p.K();
        Cp[i] = p.Cp();
        typeId[i] = p.typeId();
        
        newParcel[i] = p.newParcel();
        active[i] = p.active();

        i++;
    }

    U.write();
    UCorrect.write();
    RWF.write();
    d.write();
    rho.write();
    T.write();
    K.write();
    Cp.write();
    typeId.write();
    
    newParcel.write();
    active.write();
    
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os, 
    const solidParcel& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.UCorrect()
            << token::SPACE << p.RWF()
            << token::SPACE << p.d()
            << token::SPACE << p.rho()
            << token::SPACE << p.T()
            << token::SPACE << p.K()
            << token::SPACE << p.Cp()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << p.active();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            
            sizeof(p.U())
            +sizeof(p.UCorrect())
            +sizeof(p.RWF())
            +sizeof(p.d())
            +sizeof(p.rho())
            +sizeof(p.T())
            +sizeof(p.K())
            +sizeof(p.Cp())
            +sizeof(p.typeId())
            
            +sizeof(p.newParcel())
            +sizeof(p.active())
        );
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
