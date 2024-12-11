/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

Class
    Foam::PackingModel

Group
    grpLagrangianIntermediateMPPICPackingSubModels

Description
    Base class for packing models.

SourceFiles
    PackingModel.C
    PackingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PackingModel_H
#define PackingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "solidParcelCloud.H"
#include "ParticleStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

// class ParticleStressModel;

// class CorrectionLimitingMethod;

template<class Type>
class AveragingMethod;

/*---------------------------------------------------------------------------*\
                         Class PackingModel Declaration
\*---------------------------------------------------------------------------*/

// template<class CloudType>
class PackingModel
{
protected:

    //- Protected data

        solidParcelCloud& cloud_;

        //- Particle stress model
        autoPtr<ParticleStressModel> particleStressModel_;


public:

    //- Runtime type information
    TypeName("packingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PackingModel,
        dictionary,
        (
            const dictionary& dict,
            solidParcelCloud& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        // PackingModel(solidParcelCloud& cloud);

        //- Construct from components
        PackingModel
        (
            const dictionary& dict,
            solidParcelCloud& cloud
            // const word& type
        );

        //- Construct copy
        // PackingModel(const PackingModel<CloudType>& cm);

        //- Construct and return a clone
        // virtual autoPtr<PackingModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~PackingModel();


    //- Selector
    static autoPtr<PackingModel> New
    (
        const dictionary& dict,
        solidParcelCloud& cloud
    );


    // Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            // solidParcelCloud::solidParcel& p,
            solidParcel& p,
            const scalar deltaT
        ) const = 0;

        virtual bool active() const = 0;

        virtual void cacheFields(const bool store);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "PackingModel.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
