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
    Foam::RanzMarshall

Group
    

Description
    The Ranz-Marshall correlation for heat transfer
    RanzMarshall correction for the calculation of Nusselt number in the condition
    Re<10000 and Ma<6

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshall_H
#define RanzMarshall_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/
class RanzMarshall
:
    public HeatTransferModel
{

public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Constructors

        //- Construct from dictionary
        RanzMarshall(const dictionary& dict, solidParcelCloud& cloud);

        //- Construct copy
        RanzMarshall(const RanzMarshall& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel> clone() const
        {
            return autoPtr<HeatTransferModel>
            (
                new RanzMarshall(*this)
            );
        }


    //- Destructor
    virtual ~RanzMarshall();


    // Member Functions

        // Evaluation
            virtual bool active() const;

            //- Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr,
                const scalar Ma
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
