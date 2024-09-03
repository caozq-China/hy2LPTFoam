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
    Foam::Kavanau

Group
    

Description
    The Kavanau correlation for heat transfer
    From "MULTIPHASE FLOWS with DROPLETS and PARTICLES"
    Consider rarefied effect P.106
    

\*---------------------------------------------------------------------------*/

#ifndef Kavanau_H
#define Kavanau_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Kavanau Declaration
\*---------------------------------------------------------------------------*/
class Kavanau
:
    public HeatTransferModel
{

public:

    //- Runtime type information
    TypeName("Kavanau");


    // Constructors

        //- Construct from dictionary
        Kavanau(const dictionary& dict, solidParcelCloud& cloud);

        //- Construct copy
        Kavanau(const Kavanau& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel> clone() const
        {
            return autoPtr<HeatTransferModel>
            (
                new Kavanau(*this)
            );
        }


    //- Destructor
    virtual ~Kavanau();


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
