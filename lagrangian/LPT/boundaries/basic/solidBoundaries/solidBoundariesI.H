/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline List< autoPtr<solidPatchBoundary> >& solidBoundaries::solidPatchBoundaryModels()
{
    return solidPatchBoundaryModels_;
}

inline const List< autoPtr<solidPatchBoundary> >& solidBoundaries::solidPatchBoundaryModels() const
{
    return solidPatchBoundaryModels_;
}

inline const List<word>& solidBoundaries::solidPatchBoundaryNames() const
{
    return solidPatchBoundaryNames_;
}

inline const List<label>& solidBoundaries::solidPatchBoundaryIds() const
{
    return solidPatchBoundaryIds_;
}

inline const List<label>& solidBoundaries::solidPatchToModelIds() const
{
    return solidPatchToModelId_;
}




inline List< autoPtr<solidCyclicBoundary> >& solidBoundaries::solidCyclicBoundaryModels()
{
    return solidCyclicBoundaryModels_;
}

inline const List< autoPtr<solidCyclicBoundary> >& solidBoundaries::solidCyclicBoundaryModels() const
{
    return solidCyclicBoundaryModels_;
}

inline const List<word>& solidBoundaries::solidCyclicBoundaryNames() const
{
    return solidCyclicBoundaryNames_;
}

inline const List<label>& solidBoundaries::solidCyclicBoundaryIds() const
{
    return solidCyclicBoundaryIds_;
}

inline const List<label>& solidBoundaries::solidCyclicBoundaryToModelIds() const
{
    return solidCyclicBoundaryToModelId_;
}



inline List< autoPtr<solidGeneralBoundary> >& solidBoundaries::solidGeneralBoundaryModels()
{
    return solidGeneralBoundaryModels_;
}

inline const List< autoPtr<solidGeneralBoundary> >& solidBoundaries::solidGeneralBoundaryModels() const
{
    return solidGeneralBoundaryModels_;
}

inline const List<word>& solidBoundaries::solidGeneralBoundaryNames() const
{
    return solidGeneralBoundaryNames_;
}

inline const List<label>& solidBoundaries::solidGeneralBoundaryIds() const
{
    return solidGeneralBoundaryIds_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
