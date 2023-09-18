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
#include "constants.H"
using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const Foam::IOdictionary& Foam::solidParticleCouplingCloud::mppicProperties() const
{
    return mppicPropertiesDict_;
}

inline const Foam::List<Foam::word>& Foam::solidParticleCouplingCloud::typeIdSolidList() const
{
    return typeIdSolidList_;
}

inline const Foam::List<Foam::word>& Foam::solidParticleCouplingCloud::materialList() const
{
    return materialList_;
}

inline Foam::scalar Foam::solidParticleCouplingCloud::nSolidParticles() const
{
    return nSolidParticles_;
}

inline Foam::scalar Foam::solidParticleCouplingCloud::solidWeightingAmplificationFactor() const
{
    return solidWeightFactor_;
}

inline Foam::scalar Foam::solidParticleCouplingCloud::criticalVolFrac() const
{
    return criticalVolFrac_;
}


inline Foam::dsmcCloud* Foam::solidParticleCouplingCloud::dsmcCloudReference() const
{
    return dsmcCloudReference_;
}

inline const Foam::solidBinaryCollisionModel&  Foam::solidParticleCouplingCloud::solidBinaryCollision() const
{
    return solidBinaryCollisionModel_;
}

inline Foam::solidBinaryCollisionModel& Foam::solidParticleCouplingCloud::solidBinaryCollision()
{
    return solidBinaryCollisionModel_();
}

inline const Foam::solidCollisionDetection&  Foam::solidParticleCouplingCloud::solidCollisionDetectionModel() const
{
    return solidCollisionDetectionModel_;
}
            
inline Foam::solidCollisionDetection& Foam::solidParticleCouplingCloud::solidCollisionDetectionModel()
{
    return solidCollisionDetectionModel_();
}

inline Foam::volScalarField& Foam::solidParticleCouplingCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}

inline Foam::volScalarField& Foam::solidParticleCouplingCloud::cRsigmaPPMax()
{
    return cRsigmaPPMax_;
}

inline Foam::volScalarField& Foam::solidParticleCouplingCloud::volumeFrac()
{
    return volumeFraction_;
}

// inline Foam::volScalarField& Foam::solidParticleCouplingCloud::massBelowLunarSurface()
// {
//     return massBelowLunarSurface_;
// }
// 
// inline Foam::volScalarField& Foam::solidParticleCouplingCloud::massErosionRate()
// {
//     return massErosionRate_;
// }

inline const Foam::List<Foam::solidParticleCoupling::constantProperties>& Foam::solidParticleCouplingCloud::constSolidProps() const
{
    return constSolidProps_;
}

inline const Foam::solidParticleCoupling::constantProperties& Foam::solidParticleCouplingCloud::constSolidProps
(
    const label typeIdSolid
) const
{
    if (typeIdSolid < 0 || typeIdSolid >= constSolidProps_.size())
    {
        FatalErrorInFunction
            << "solid constantProperties for requested typeId index "
            << typeIdSolid << " do not exist" << nl
            << abort(FatalError);
    }

    return constSolidProps_[typeIdSolid];
}

inline Foam::Random& Foam::solidParticleCouplingCloud::rndGenS()
{
    return rndGenS_;
}


inline bool Foam::solidParticleCouplingCloud::enableMPPICMethod()
{
    return enableMPPICMethod_;
}

inline const  Foam::solidBoundaries&  Foam::solidParticleCouplingCloud::boundaries() const
{
    return  boundaries_;
}

inline  Foam::solidBoundaries&  Foam::solidParticleCouplingCloud::boundaries()
{
    return  boundaries_;
}

inline Foam::solidBoundaryMeasurements& Foam::solidParticleCouplingCloud::boundaryFluxMeasurements()
{
    return solidBoundaryMeas_;
}

inline const Foam::solidBoundaryMeasurements& Foam::solidParticleCouplingCloud::boundaryFluxMeasurements() const
{
    return solidBoundaryMeas_;
}

inline Foam::solidCellMeasurements& Foam::solidParticleCouplingCloud::solidCellPropMeasurements()
{
    return solidCellMeas_;
}


inline const Foam::solidCellMeasurements& Foam::solidParticleCouplingCloud::solidCellPropMeasurements() const
{
    return solidCellMeas_;
}


inline const Foam::solidControllers& Foam::solidParticleCouplingCloud::controllers() const
{
    return controllers_;
}

inline Foam::solidControllers& Foam::solidParticleCouplingCloud::controllers()
{
    return controllers_;
}

inline const Foam::InterphaseCoupling& Foam::solidParticleCouplingCloud::interphaseCoupling() const
{
    return interphaseCouplingModel_;
}

inline Foam::InterphaseCoupling& Foam::solidParticleCouplingCloud::interphaseCoupling()
{
    return interphaseCouplingModel_();
}


inline const Foam::solidPhaseChangeModel& Foam::solidParticleCouplingCloud::solidPhaseChange() const
{
    return solidPhaseChange_;
}

inline Foam::solidPhaseChangeModel& Foam::solidParticleCouplingCloud::solidPhaseChange()
{
    return solidPhaseChange_();
}

// inline const Foam::surfaceErosionModel& Foam::solidParticleCouplingCloud::surfaceErosionModels() const
// {
//     return surfaceErosionModel_;
// }
// 
// inline Foam::surfaceErosionModel& Foam::solidParticleCouplingCloud::surfaceErosionModels()
// {
//     return surfaceErosionModel_();
// }

inline const Foam::DynamicList<Foam::label>& Foam::solidParticleCouplingCloud::cellSolidPackingOccupancyIds() const
{
    return cellSolidPackingOccupIds_;
}

inline const Foam::DynamicList<Foam::label>& Foam::solidParticleCouplingCloud::cellSolidCollisionOccupancyIds() const
{
    return cellSolidCollisionOccupIds_;
}

inline const Foam::mppicDampingModel& Foam::solidParticleCouplingCloud::mppicDampingModels() const
{
    return mppicDampingModel_;
}

inline Foam::mppicDampingModel& Foam::solidParticleCouplingCloud::mppicDampingModels()
{
    return mppicDampingModel_();
}

inline const Foam::mppicPackingModel& Foam::solidParticleCouplingCloud::mppicPackingModels() const
{
    return mppicPackingModel_;
}

inline Foam::mppicPackingModel& Foam::solidParticleCouplingCloud::mppicPackingModels()
{
    return mppicPackingModel_();
}

inline const Foam::mppicIsotropyModel& Foam::solidParticleCouplingCloud::mppicIsotropyModels() const
{
    return mppicIsotropyModel_;
}

inline Foam::mppicIsotropyModel& Foam::solidParticleCouplingCloud::mppicIsotropyModels()
{
    return mppicIsotropyModel_();
}

inline void Foam::solidParticleCouplingCloud::clearSolidParticles()
{
    return IDLList<solidParticleCoupling>::clear();
}

// ************************************************************************* //
