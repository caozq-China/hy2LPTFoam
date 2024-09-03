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
    Foam::DampingModel

Group
    grpLagrangianIntermediateMPPICDampingSubModels

Description
    Base class for collisional damping models.

SourceFiles
    DampingModel.C
    DampingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DampingModel_H
#define DampingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
// #include "CloudSubModelBase.H"
#include "solidParcelCloud.H"
#include "TimeScaleModel.H"
#include "AveragingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class TimeScaleModel;

/*---------------------------------------------------------------------------*\
                         Class DampingModel Declaration
\*---------------------------------------------------------------------------*/

class DampingModel
{
protected:

    // Protected data
        solidParcelCloud&  cloud_;

        //- Time scale model
        autoPtr<TimeScaleModel> timeScaleModel_;


public:

    //- Runtime type information
    TypeName("dampingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DampingModel,
        dictionary,
        (
            const dictionary& dict,
            solidParcelCloud& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        // DampingModel(solidParcelCloud& cloud);

        //- Construct from components
        DampingModel
        (
            const dictionary& dict,
            solidParcelCloud& cloud
            // const word& type
        );

        //- Construct copy
        // DampingModel(const DampingModel& cm);

        //- Construct and return a clone
        // virtual autoPtr<DampingModel> clone() const = 0;


    //- Destructor
    virtual ~DampingModel();


    //- Selector
    static autoPtr<DampingModel> New
    (
        const dictionary& dict,
        solidParcelCloud& cloud
    );


    // Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            solidParcel& p,
            const scalar deltaT
        ) const = 0;

        virtual bool active() const = 0;

        // virtual void calUcorrect();

        virtual void cacheFields(const bool store) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "DampingModel.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
