/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::HeatTransferModel

Group

Description
    Templated heat transfer model class

SourceFiles
    HeatTransferModel.C
    HeatTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferModel_H
#define HeatTransferModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class solidParcelCloud;
/*---------------------------------------------------------------------------*\
                       Class HeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class HeatTransferModel
{
protected:
    // Protected data

        //- Apply Bird's correction to the htc
        const Switch BirdCorrection_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HeatTransferModel,
        dictionary,
        (
            const dictionary& dict,
            solidParcelCloud& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        HeatTransferModel(solidParcelCloud& cloud);

        //- Construct from dictionary
        HeatTransferModel
        (
            const dictionary& dict,
            solidParcelCloud& cloud,
            const word& type
        );

        //- Construct copy
        HeatTransferModel(const HeatTransferModel& htm);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel> clone() const = 0;


    //- Destructor
    virtual ~HeatTransferModel();


    //- Selector
    static autoPtr<HeatTransferModel> New
    (
        const dictionary& dict,
        solidParcelCloud& cloud
    );


    // Member Functions

        // Flag

            virtual bool active() const = 0;

        // Access

            //- Return the Bird htc correction flag
            const Switch& BirdCorrection() const;


        // Evaluation

            //- Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr,
                const scalar Ma
            ) const = 0;

            //- Return heat transfer coefficient
            virtual scalar htc
            (
                const scalar dp,
                const scalar Re,
                const scalar Pr,
                const scalar kappa,
                const scalar NCpW,
                const scalar Ma
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
