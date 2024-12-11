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


Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- Access
inline const List< List < vectorField > >& boundaryMeasurements::momentumBF() const
{
    return momentumBF_;
}

inline const List< List < vectorField > >& boundaryMeasurements::UMeanBF() const
{
    return UMeanBF_;
}

inline const List< List < vectorField > >& boundaryMeasurements::fDBF() const
{
    return fDBF_;
}

inline const List< List < scalarField > >& boundaryMeasurements::rhoNBF() const
{
    return rhoNBF_;
}

inline const List< List < scalarField > >& boundaryMeasurements::rhoMBF() const
{
    return rhoMBF_;
}

inline const List< List < scalarField > >& boundaryMeasurements::qBF() const
{
    return qBF_;
}


// Edit
inline List< List < vectorField > >& boundaryMeasurements::momentumBF()
{
    return momentumBF_;
}

inline List< List < vectorField > >& boundaryMeasurements::UMeanBF()
{
    return UMeanBF_;
}

inline List< List < vectorField > >& boundaryMeasurements::fDBF()
{
    return fDBF_;
}

inline List< List < scalarField > >& boundaryMeasurements::rhoNBF()
{
    return rhoNBF_;
}

inline List< List < scalarField > >& boundaryMeasurements::rhoMBF()
{
    return rhoMBF_;
}

inline List< List < scalarField > >& boundaryMeasurements::qBF()
{
    return qBF_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
