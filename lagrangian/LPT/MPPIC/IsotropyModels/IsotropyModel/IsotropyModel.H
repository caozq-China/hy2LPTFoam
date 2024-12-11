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
    Foam::IsotropyModel

Group
    grpLagrangianIntermediateMPPICIsotropySubModels

Description
    Base class for collisional return-to-isotropy models.

SourceFiles
    IsotropyModel.C
    IsotropyModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef IsotropyModel_H
#define IsotropyModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "solidParcelCloud.H"
#include "TimeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

// class TimeScaleModel;

/*---------------------------------------------------------------------------*\
                        Class IsotropyModel Declaration
\*---------------------------------------------------------------------------*/

class IsotropyModel
{
protected:

        solidParcelCloud& cloud_;

        //- Time scale model
        autoPtr<TimeScaleModel> timeScaleModel_;


public:

    //- Runtime type information
    TypeName("isotropyModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        IsotropyModel,
        dictionary,
        (
            const dictionary& dict,
            solidParcelCloud& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        // IsotropyModel(solidParcelCloud& cloud);

        //- Construct from components
        IsotropyModel
        (
            const dictionary& dict,
            solidParcelCloud& cloud
            // const word& type
        );

        //- Construct a copy
        // IsotropyModel(const IsotropyModel<CloudType>& cm);

        //- Construct and return a clone
        // virtual autoPtr<IsotropyModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~IsotropyModel();


    //- Selector
    static autoPtr<IsotropyModel> New
    (
        const dictionary& dict,
        solidParcelCloud& cloud
    );


    //- Member Functions

        //- Calculate velocities
        virtual void calculate() = 0;

        virtual bool active() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "IsotropyModel.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
