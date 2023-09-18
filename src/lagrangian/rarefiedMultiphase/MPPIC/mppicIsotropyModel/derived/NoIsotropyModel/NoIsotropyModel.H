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
    Foam::IsotropyModels::NoIsotropy

Description

SourceFiles
    NoIsotropy.C

\*---------------------------------------------------------------------------*/

#ifndef NoIsotropyModel_H
#define NoIsotropyModel_H

#include "mppicIsotropyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace mppicIsotropyModel
// {

/*---------------------------------------------------------------------------*\
                         Class NoIsotropy Declaration
\*---------------------------------------------------------------------------*/

// template<class CloudType>
class NoIsotropyModel
:
    public mppicIsotropyModel
{
public:

    //- Runtime type information
    TypeName("NoIsotropyModel");

    // Constructors

        //- Construct from components
        NoIsotropyModel
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
        );

        //- Construct copy
//         NoIsotropy(const NoIsotropy<CloudType>& cm);

        //- Construct and return a clone
//         virtual autoPtr<IsotropyModel<CloudType> > clone() const
//         {
//             return autoPtr<IsotropyModel<CloudType> >
//             (
//                 new NoIsotropy<CloudType>(*this)
//             );
//         }


    //- Destructor
    virtual ~NoIsotropyModel();


    //- Member Functions

        //- calculate velocities
        virtual void calculate();

        //- Return the model 'active' status
        virtual bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace IsotropyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
