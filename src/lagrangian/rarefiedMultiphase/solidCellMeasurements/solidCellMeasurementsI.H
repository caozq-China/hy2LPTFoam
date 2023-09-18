/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// inline const scalarField& cellMeasurements::collisionSeparation() const
// {
//     return collisionSeparation_;
// }
// 
// 
// inline scalarField& cellMeasurements::collisionSeparation()
// {
//     return collisionSeparation_;
// }


inline const scalarField& solidCellMeasurements::nColls() const
{
    return nColls_;
}


inline scalarField& solidCellMeasurements::nColls()
{
    return nColls_;
}

inline const scalarField& solidCellMeasurements::nCollsTotal() const
{
    return nCollsTotal_;
}


inline scalarField& solidCellMeasurements::nCollsTotal()
{
    return nCollsTotal_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
