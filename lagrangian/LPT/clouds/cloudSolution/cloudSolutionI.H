/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

inline const Foam::fvMesh& Foam::cloudSolution::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::cloudSolution::dict() const
{
    return dict_;
}


inline const Foam::Switch Foam::cloudSolution::active() const
{
    return active_;
}

inline const Foam::Switch Foam::cloudSolution::smoothTheta() const
{
    return smoothTheta_;
}

// inline Foam::Switch& Foam::cloudSolution::smoothTheta()
// {
//     return smoothTheta_;
// }

inline const Foam::Switch Foam::cloudSolution::smoothSourceTerms() const
{
    return smoothSourceTerms_;
}

inline const Foam::Switch Foam::cloudSolution::useGasFilter() const
{
    return useGasFilter_;
}

inline const Foam::Switch Foam::cloudSolution::useGasFilterImplicit() const
{
    return useGasFilterImplicit_;
}


inline const Foam::dictionary& Foam::cloudSolution::sourceTermDict() const
{
    return dict_.subDict("sourceTerms");
}


inline const Foam::dictionary& Foam::cloudSolution::interpolationSchemes() const
{
    return dict_.subDict("interpolationSchemes");
}


inline const Foam::dictionary& Foam::cloudSolution::integrationSchemes() const
{
    return dict_.subDict("integrationSchemes");
}


inline Foam::scalar Foam::cloudSolution::maxCo() const
{
    return maxCo_;
}


inline Foam::label Foam::cloudSolution::iter() const
{
    return iter_;
}


inline Foam::label Foam::cloudSolution::nextIter()
{
    return ++iter_;
}


inline Foam::scalar Foam::cloudSolution::trackTime() const
{
    return trackTime_;
}


inline Foam::scalar Foam::cloudSolution::deltaTMax() const
{
    return deltaTMax_;
}


inline Foam::Switch& Foam::cloudSolution::coupled()
{
    return coupled_;
}


inline const Foam::Switch Foam::cloudSolution::coupled() const
{
    return coupled_;
}


inline const Foam::Switch Foam::cloudSolution::cellValueSourceCorrection() const
{
    return cellValueSourceCorrection_;
}


// ************************************************************************* //
