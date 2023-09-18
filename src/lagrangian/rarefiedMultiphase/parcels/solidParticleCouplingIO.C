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

#include "solidParticleCoupling.H"
#include "solidParticleCouplingCloud.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticleCoupling::solidParticleCoupling
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    U_(Zero),
    omega_(Zero),
    UCorrect_(Zero),
    F_(Zero),
    D_(0.0),
    T_(0.0),
    RWF_(1.0),
    CzRatio_(0.0),
    numSeq_(-1),
    typeID_(-1),
    phaseState_(0),
    newSolidParticle_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            is >> omega_;
            is >> UCorrect_;
            is >> F_;
            D_ = readScalar(is);
            T_ = readScalar(is);
            RWF_ = readScalar(is);
            CzRatio_ = readScalar(is);
            numSeq_ = readLabel(is);
            typeID_= readLabel(is);
            phaseState_= readLabel(is);
            newSolidParticle_= readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
                + sizeof(omega_)
                + sizeof(UCorrect_)
                + sizeof(F_)
                + sizeof(D_)
                + sizeof(T_)
                + sizeof(RWF_)
                + sizeof(CzRatio_)
                + sizeof(numSeq_)
                + sizeof(typeID_)
                + sizeof(phaseState_)
                + sizeof(newSolidParticle_)
            );
        }
    }

     // Check state of Istream
    is.check(FUNCTION_NAME);
}


void Foam::solidParticleCoupling::readFields(Cloud<solidParticleCoupling>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);
    
    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<vector> omega(c.fieldIOobject("omega", IOobject::MUST_READ));
    c.checkFieldIOobject(c, omega);
    
    IOField<vector> UCorrect(c.fieldIOobject("UCorrect", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UCorrect);
    
    IOField<vector> F(c.fieldIOobject("F", IOobject::MUST_READ));
    c.checkFieldIOobject(c, F);
    
    IOField<scalar> D(c.fieldIOobject("D", IOobject::MUST_READ));
    c.checkFieldIOobject(c, D);
    
    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ));
    c.checkFieldIOobject(c, T);
    
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::MUST_READ));
    c.checkFieldIOobject(c, RWF);
    
    IOField<scalar> CzRatio(c.fieldIOobject("crystallizationRatio", IOobject::MUST_READ));
    c.checkFieldIOobject(c, CzRatio);

    IOField<label> numSeq(c.fieldIOobject("numSeq", IOobject::MUST_READ));
    c.checkFieldIOobject(c, numSeq);
    
    IOField<label> typeID(c.fieldIOobject("typeID", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeID);
    
    IOField<label> phaseState(c.fieldIOobject("phaseState", IOobject::MUST_READ));
    c.checkFieldIOobject(c, phaseState);
    
    IOField<label> newSolidParticle(c.fieldIOobject("newSolidParticle", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newSolidParticle);
    
    label i = 0;
    forAllIter(solidParticleCouplingCloud, c, iter)
    {
        solidParticleCoupling& p = iter();

        p.U_ = U[i];
        p.omega_ = omega[i];
        p.UCorrect_ = UCorrect[i];
        p.F_ = F[i];
        p.D_ = D[i];
        p.T_ = T[i];
        p.RWF_ = RWF[i];
        p.CzRatio_ = CzRatio[i];
        p.numSeq_ = numSeq[i];
        p.typeID_ = typeID[i];
        p.phaseState_ = phaseState[i];
        p.newSolidParticle_ = newSolidParticle[i];
        i++;
    }
}


void Foam::solidParticleCoupling::writeFields(const Cloud<solidParticleCoupling>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> omega(c.fieldIOobject("omega", IOobject::NO_READ), np);
    IOField<vector> UCorrect(c.fieldIOobject("UCorrect", IOobject::NO_READ), np);
    IOField<vector> F(c.fieldIOobject("F", IOobject::NO_READ), np);
    IOField<scalar> D(c.fieldIOobject("D", IOobject::NO_READ), np);
    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<scalar> CzRatio(c.fieldIOobject("crystallizationRatio", IOobject::NO_READ), np);
    IOField<label> numSeq(c.fieldIOobject("numSeq", IOobject::NO_READ), np);
    IOField<label> typeID(c.fieldIOobject("typeID", IOobject::NO_READ), np);
    IOField<label> phaseState(c.fieldIOobject("phaseState", IOobject::NO_READ), np);
    IOField<label> newSolidParticle(c.fieldIOobject("newSolidParticle", IOobject::NO_READ), np);


    label i = 0;
    forAllConstIter(solidParticleCouplingCloud, c, iter)
    {
        const solidParticleCoupling& p = iter();
        
        U[i] = p.U();
        omega[i] = p.omega();
        UCorrect[i] = p.UCorrect();
        F[i] = p.F();
        D[i] = p.D();
        T[i] = p.T();
        RWF[i] = p.RWF();
        CzRatio[i] = p.CzRatio();
        numSeq[i] = p.numSeq();
        typeID[i] = p.typeID();
        phaseState[i] = p.phaseState();
        newSolidParticle[i] = p.newSolidParticle();
        i++;
    }

    U.write();
    omega.write();
    UCorrect.write();
    F.write();
    D.write();
    T.write();
    RWF.write();
    CzRatio.write();
    numSeq.write();
    typeID.write();
    phaseState.write();
    newSolidParticle.write();

}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os, 
    const solidParticleCoupling& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)   
            << token::SPACE << p.U()
            << token::SPACE << p.omega()
            << token::SPACE << p.UCorrect()
            << token::SPACE << p.F()
            << token::SPACE << p.D()
            << token::SPACE << p.T()
            << token::SPACE << p.RWF()
            << token::SPACE << p.CzRatio()
            << token::SPACE << p.numSeq()
            << token::SPACE << p.typeID()
            << token::SPACE << p.phaseState()
            << token::SPACE << p.newSolidParticle();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            
            sizeof(p.U())
            +sizeof(p.omega())
            + sizeof(p.UCorrect())
            + sizeof(p.F())
            + sizeof(p.D())
            + sizeof(p.T())
            + sizeof(p.RWF())
            + sizeof(p.CzRatio())
            + sizeof(p.numSeq())
            + sizeof(p.typeID()) 
            + sizeof(p.phaseState())
            + sizeof(p.newSolidParticle())
        );
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
