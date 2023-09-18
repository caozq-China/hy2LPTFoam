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

#include "xParticle.H"
#include "xParticleCloud.H"
#include "IOstreams.H"
// #include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// const std::size_t Foam::xParticle::sizeofFields
// (
//     sizeof(xParticle) - offsetof(xParticle,typeID_)
// );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::xParticle::xParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    typeID_(0),
    newParticle_(0),
    T_(0),
    RWF_(0),
    U_(Zero),
    // omega_(Zero),
    f_(Zero),
    angularMomentum_(Zero),
    torque_(Zero),
    collisionRecords_()
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> typeID_
                >> newParticle_
                >> T_
                >> RWF_
                >> U_
                // >> omega_
                >> f_
                >> angularMomentum_
                >> torque_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size
            is.beginRawRead();
            readRawLabel(is, &typeID_);
            readRawLabel(is, &newParticle_);
            readRawScalar(is, &T_);
            readRawScalar(is, &RWF_);
            readRawScalar(is, U_.data(), vector::nComponents);
            // readRawScalar(is, omega_.data(), vector::nComponents);
            readRawScalar(is, f_.data(), vector::nComponents);
            readRawScalar(is, angularMomentum_.data(), vector::nComponents);
            readRawScalar(is, torque_.data(), vector::nComponents);
            is.endRawRead();

        }
        else
        {
            is.read(reinterpret_cast<char*>(&U_),
//                     sizeofFields);
                    sizeof(typeID_)
                   +sizeof(newParticle_)
                   +sizeof(T_)
                   +sizeof(RWF_)
                   +sizeof(U_)
                   // +sizeof(omega_)
                   +sizeof(f_)
                   +sizeof(angularMomentum_)
                   +sizeof(torque_)
               );
        }
        is >> collisionRecords_;
    }

     // Check state of Istream
    is.check(FUNCTION_NAME);
}


void Foam::xParticle::readFields(Cloud<xParticle>& c)
{
    const bool readOnProc = c.size();

    particle::readFields(c);

    IOField<label> typeID(c.fieldIOobject("typeID", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, typeID);

    IOField<label> newParticle(c.fieldIOobject("newParticle", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, newParticle);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, T);

    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, RWF);
    
    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, U);

    // IOField<vector> omega(c.fieldIOobject("omega", IOobject::MUST_READ),readOnProc);
    // c.checkFieldIOobject(c, omega);

    IOField<vector> f(c.fieldIOobject("f", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, f);

    IOField<vector> angularMomentum(c.fieldIOobject("angularMomentum", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, angularMomentum);

    IOField<vector> torque(c.fieldIOobject("torque", IOobject::MUST_READ),readOnProc);
    c.checkFieldIOobject(c, torque);

    //- collsion records
    labelFieldCompactIOField collisionRecordsPairAccessed
    (
        c.fieldIOobject("collisionRecordsPairAccessed", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairAccessed);

    labelFieldCompactIOField collisionRecordsPairOrigProcOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigProcOfOther",
            IOobject::MUST_READ
        ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairOrigProcOfOther);

    labelFieldCompactIOField collisionRecordsPairOrigIdOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigIdOfOther",
            IOobject::MUST_READ
        ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairOrigProcOfOther);

    pairDataFieldCompactIOField collisionRecordsPairData
    (
        c.fieldIOobject("collisionRecordsPairData", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairData);

    labelFieldCompactIOField collisionRecordsWallAccessed
    (
        c.fieldIOobject("collisionRecordsWallAccessed", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallAccessed);

    vectorFieldCompactIOField collisionRecordsWallPRel
    (
        c.fieldIOobject("collisionRecordsWallPRel", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallPRel);

    wallDataFieldCompactIOField collisionRecordsWallData
    (
        c.fieldIOobject("collisionRecordsWallData", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallData);

    label i = 0;
    for(xParticle& p : c)
    {
        p.typeID_ = typeID[i];
        p.newParticle_ = newParticle[i];
        p.T_ = T[i];
        p.RWF_ = RWF[i];
        p.U_ = U[i];
        // p.omega_ = omega[i];
        p.f_ = f[i];
        p.angularMomentum_ = angularMomentum[i];
        p.torque_ = torque[i];

        p.collisionRecords_ = collisionRecordList
        (
            collisionRecordsPairAccessed[i],
            collisionRecordsPairOrigProcOfOther[i],
            collisionRecordsPairOrigIdOfOther[i],
            collisionRecordsPairData[i],
            collisionRecordsWallAccessed[i],
            collisionRecordsWallPRel[i],
            collisionRecordsWallData[i]
        );
        ++i;
    }
}


void Foam::xParticle::writeFields(const Cloud<xParticle>& c)
{
    particle::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<label> typeID(c.fieldIOobject("typeID", IOobject::NO_READ), np);
    IOField<label> newParticle(c.fieldIOobject("newParticle", IOobject::NO_READ), np);
    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    // IOField<vector> omega(c.fieldIOobject("omega", IOobject::NO_READ), np);
    IOField<vector> f(c.fieldIOobject("f", IOobject::NO_READ), np);
    IOField<vector> angularMomentum(c.fieldIOobject("angularMomentum", IOobject::NO_READ), np);
    IOField<vector> torque(c.fieldIOobject("torque", IOobject::NO_READ), np);

    //- collision records
    labelFieldCompactIOField collisionRecordsPairAccessed
    (
        c.fieldIOobject("collisionRecordsPairAccessed", IOobject::NO_READ),
        np
    );
    labelFieldCompactIOField collisionRecordsPairOrigProcOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigProcOfOther",
            IOobject::NO_READ
        ),
        np
    );
    labelFieldCompactIOField collisionRecordsPairOrigIdOfOther
    (
        c.fieldIOobject("collisionRecordsPairOrigIdOfOther", IOobject::NO_READ),
        np
    );
    pairDataFieldCompactIOField collisionRecordsPairData
    (
        c.fieldIOobject("collisionRecordsPairData", IOobject::NO_READ),
        np
    );
    labelFieldCompactIOField collisionRecordsWallAccessed
    (
        c.fieldIOobject("collisionRecordsWallAccessed", IOobject::NO_READ),
        np
    );
    vectorFieldCompactIOField collisionRecordsWallPRel
    (
        c.fieldIOobject("collisionRecordsWallPRel", IOobject::NO_READ),
        np
    );
    wallDataFieldCompactIOField collisionRecordsWallData
    (
        c.fieldIOobject("collisionRecordsWallData", IOobject::NO_READ),
        np
    );


    label i = 0;
    for(const xParticle& p : c)
    {
        typeID[i] = p.typeID_;
        newParticle[i] = p.newParticle_;
        T[i] = p.T_;
        RWF[i] = p.RWF_;
        U[i] = p.U_;
        // omega[i] = p.omega_;
        f[i] = p.f_;
        angularMomentum[i] = p.angularMomentum_;
        torque[i] = p.torque_;
        collisionRecordsPairAccessed[i] = p.collisionRecords().pairAccessed();
        collisionRecordsPairOrigProcOfOther[i] = p.collisionRecords().pairOrigProcOfOther();
        collisionRecordsPairOrigIdOfOther[i] = p.collisionRecords().pairOrigIdOfOther();
        collisionRecordsPairData[i] = p.collisionRecords().pairData();
        collisionRecordsWallAccessed[i] = p.collisionRecords().wallAccessed();
        collisionRecordsWallPRel[i] = p.collisionRecords().wallPRel();
        collisionRecordsWallData[i] = p.collisionRecords().wallData();
        ++i;
    }

    typeID.write(writeOnProc);
    newParticle.write(writeOnProc);
    T.write(writeOnProc);
    RWF.write(writeOnProc);
    U.write(writeOnProc);
    // omega.write(writeOnProc);
    f.write(writeOnProc);
    angularMomentum.write(writeOnProc);
    torque.write(writeOnProc);

    collisionRecordsPairAccessed.write(writeOnProc);
    collisionRecordsPairOrigProcOfOther.write(writeOnProc);
    collisionRecordsPairOrigIdOfOther.write(writeOnProc);
    collisionRecordsPairData.write(writeOnProc);
    collisionRecordsWallAccessed.write(writeOnProc);
    collisionRecordsWallPRel.write(writeOnProc);
    collisionRecordsWallData.write(writeOnProc);

}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os, 
    const xParticle& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.typeID_
            << token::SPACE << p.newParticle_
            << token::SPACE << p.T_
            << token::SPACE << p.RWF_
            << token::SPACE << p.U_
            // << token::SPACE << p.omega_
            << token::SPACE << p.f_
            << token::SPACE << p.angularMomentum_
            << token::SPACE << p.torque_
            << token::SPACE << p.collisionRecords_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            sizeof(p.typeID_)
           +sizeof(p.newParticle_)
           +sizeof(p.T_)
           +sizeof(p.RWF_)
           +sizeof(p.U_)
           // +sizeof(p.omega_)
           +sizeof(p.f_)
           +sizeof(p.angularMomentum_)
           +sizeof(p.torque_)

        );
        os << p.collisionRecords();
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
