/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Description
    Base class for collisional damping models.

SourceFiles
    DampingModel.C
    DampingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef mppicDampingModel_H
#define mppicDampingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "solidParticleCouplingCloud.H"
#include "mppicTimeScaleModels.H"
// #include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class mppicTimeScaleModels;

/*---------------------------------------------------------------------------*\
                         Class mppicDampingModel Declaration
\*---------------------------------------------------------------------------*/
class mppicDampingModel
// :
//     public CloudSubModelBase<CloudType>
{
protected:

    // Protected data
    
//         word regionName_;
        
//         label regionID_;
        
        solidParticleCouplingCloud& spc_;

        //- Time scale model
        autoPtr<mppicTimeScaleModels> mppicTimeScaleModel_;


public:

    //- Runtime type information
    TypeName("mppicDampingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mppicDampingModel,
        dictionary,
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
//             CloudType& owner
        ),
        (spc,dict)
    );


    // Constructors

        //- Construct null from owner
//         mppicDampingModel(CloudType& owner);

        //- Construct from components
        mppicDampingModel
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
            
//             const word& type
        );

        //- Construct copy
//         mppicDampingModel(const mppicDampingModel<CloudType>& cm);

        //- Construct and return a clone
//         virtual autoPtr<mppicDampingModel<CloudType> > clone() const = 0;


    //- Destructor
    virtual ~mppicDampingModel();


    //- Selector
//     static autoPtr<mppicDampingModel<CloudType> > New
//     (
//         const dictionary& dict,
//         CloudType& owner
//     );
    static autoPtr<mppicDampingModel> New
    (
        solidParticleCouplingCloud& spc,
        const dictionary& dict
//         CloudType& owner
    );


    // Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            solidParticleCoupling& p,
            const scalar deltaT
        ) const = 0;
        
        virtual bool active() const = 0;
        
//         virtual const word& regionName() const;
        
        //- Return the control zone cells
//         virtual const labelList& controlZone() const;
        
        virtual void calUcorrect();
        
        //- Calculate the damping time scales
        virtual void cacheFields(const bool store) = 0;;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
