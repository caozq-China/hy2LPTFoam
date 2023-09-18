/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

Description
    Base class for collisional return-to-isotropy models.

SourceFiles
    IsotropyModel.C
    IsotropyModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef mppicIsotropyModel_H
#define mppicIsotropyModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "solidParticleCouplingCloud.H"
#include "mppicTimeScaleModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class mppicTimeScaleModels;

/*---------------------------------------------------------------------------*\
                        Class mppicIsotropyModels Declaration
\*---------------------------------------------------------------------------*/

// template<class CloudType>
class mppicIsotropyModel
// :
//     public CloudSubModelBase<CloudType>
{
protected:
    
//         word regionName_;
        
//         label regionID_;
        
        solidParticleCouplingCloud& spc_;

        //- Time scale model
        autoPtr<mppicTimeScaleModels> mppicTimeScaleModel_;


public:

    //- Runtime type information
    TypeName("mppicIsotropyModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mppicIsotropyModel,
        dictionary,
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        ),
        (spc, dict)
    );


    // Constructors

        //- Construct null from owner
//         IsotropyModel
//         (
//             solidParticleCouplingCloud& spc
//         );

        //- Construct from components
        mppicIsotropyModel
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );

        //- Construct a copy
//         IsotropyModel(const IsotropyModel<CloudType>& cm);

        //- Construct and return a clone
//         virtual autoPtr<IsotropyModel<CloudType> > clone() const = 0;


    //- Destructor
    virtual ~mppicIsotropyModel();


    //- Selector
    static autoPtr<mppicIsotropyModel> New
    (
        solidParticleCouplingCloud& spc,
        const dictionary& dict
    );


    //- Member Functions
//         virtual const word& regionName() const;
        
        //- Return the control zone cells
//         virtual const labelList& controlZone() const;

        //- Calculate velocities
        virtual void calculate() = 0;
        
        virtual bool active() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
