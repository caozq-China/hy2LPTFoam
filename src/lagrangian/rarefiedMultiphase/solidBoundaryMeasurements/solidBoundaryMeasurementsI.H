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

// inline const List< List < scalarField > >&  solidBoundaryMeasurements::rhoNIntBF() const
// {
//     return rhoNIntBF_;
// }
// 
// inline const List< List < scalarField > >&  solidBoundaryMeasurements::rhoNElecBF() const
// {
//     return rhoNElecBF_;
// }

inline const List< List < scalarField > >& solidBoundaryMeasurements::rhoNBF() const
{
    return rhoNBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::rhoMBF() const
{
    return rhoMBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::linearKEBF() const
{
    return linearKEBF_;
}
inline const List< List < scalarField > >& solidBoundaryMeasurements::volumeFractionBF() const
{
    return volumeFractionBF_;
}
/*
inline const List< List < scalarField > >& solidBoundaryMeasurements::cXmXTBF() const
{
    return cXmXTBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::cXmBF() const
{
    return cXmBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::TXnParticleBF() const
{
    return TXnParticleBF_;
}
*/
// inline const List< List < scalarField > >& solidBoundaryMeasurements::mccSpeciesBF() const
// {
//     return mccSpeciesBF_;
// }

inline const List< List < vectorField > >& solidBoundaryMeasurements::momentumBF() const
{
    return momentumBF_;
}

inline const List< List < vectorField > >& solidBoundaryMeasurements::UMeanBF() const
{
    return UMeanBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::numberFluxBF() const
{
    return numberFluxBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::massFluxBF() const
{
    return massFluxBF_;
}

inline const List< List < scalarField > >& solidBoundaryMeasurements::conductiveHeatFluxBF() const
{
    return conductiveHeatFluxBF_;
}
// inline const List< List < scalarField > >& solidBoundaryMeasurements::rotationalEBF() const
// {
//     return rotationalEBF_;
// }
// 
// inline const List< List < scalarField > >& solidBoundaryMeasurements::rotationalDofBF() const
// {
//     return rotationalDofBF_;
// }
// 
// inline const List< List < scalarField > >& solidBoundaryMeasurements::vibrationalEBF() const
// {
//     return vibrationalEBF_;
// }
// 
// inline const List< List < scalarField > >& solidBoundaryMeasurements::electronicEBF() const
// {
//     return electronicEBF_;
// }
// 
// inline const List< List < scalarField > >& solidBoundaryMeasurements::qBF() const
// {
//     return qBF_;
// }
// 
// inline const List< List < vectorField > >& solidBoundaryMeasurements::fDBF() const
// {
//     return fDBF_;
// }

// - Edit

// inline List< List < scalarField > >& solidBoundaryMeasurements::rhoNIntBF()
// {
//     return rhoNIntBF_;
// }
// 
// inline List< List < scalarField > >& solidBoundaryMeasurements::rhoNElecBF()
// {
//     return rhoNElecBF_;
// }

inline List< List < scalarField > >& solidBoundaryMeasurements::rhoNBF()
{
    return rhoNBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::rhoMBF()
{
    return rhoMBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::linearKEBF()
{
    return linearKEBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::volumeFractionBF()
{
    return volumeFractionBF_;
}
/*
inline List< List < scalarField > >& solidBoundaryMeasurements::cXmXTBF()
{
    return cXmXTBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::cXmBF()
{
    return cXmBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::TXnParticleBF()
{
    return TXnParticleBF_;
}
*/
// inline List< List < scalarField > >& solidBoundaryMeasurements::mccSpeciesBF()
// {
//     return mccSpeciesBF_;
// }

inline List< List < vectorField > >& solidBoundaryMeasurements::momentumBF()
{
    return momentumBF_;
}

inline List< List < vectorField > >& solidBoundaryMeasurements::UMeanBF()
{
    return UMeanBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::numberFluxBF()
{
    return numberFluxBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::massFluxBF()
{
    return massFluxBF_;
}

inline List< List < scalarField > >& solidBoundaryMeasurements::conductiveHeatFluxBF()
{
    return conductiveHeatFluxBF_;
}
// inline List< List < scalarField > >& solidBoundaryMeasurements::rotationalEBF()
// {
//     return rotationalEBF_;
// }
// 
// inline List< List < scalarField > >& solidBoundaryMeasurements::rotationalDofBF()
// {
//     return rotationalDofBF_;
// }
// 
// inline List< List < scalarField > >& solidBoundaryMeasurements::vibrationalEBF()
// {
//     return vibrationalEBF_;
// }
// 
// inline List< List < scalarField > >& solidBoundaryMeasurements::electronicEBF()
// {
//     return electronicEBF_;
// }
// 
// inline List< List < scalarField > >& solidBoundaryMeasurements::qBF()
// {
//     return qBF_;
// }
// 
// inline List< List < vectorField > >& solidBoundaryMeasurements::fDBF()
// {
//     return fDBF_;
// }


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
