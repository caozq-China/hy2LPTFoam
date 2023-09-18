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
    Foam::DampingModels::NoDamping

Description

SourceFiles
    NoDamping.C

\*---------------------------------------------------------------------------*/

#ifndef NoMppicDamping_H
#define NoMppicDamping_H

#include "mppicDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace mppicDampingModels
// {

/*---------------------------------------------------------------------------*\
                         Class NoDamping Declaration
\*---------------------------------------------------------------------------*/

class NoMppicDamping
:
    public mppicDampingModel
{
public:

    //- Runtime type information
    TypeName("NoMppicDamping");

    // Constructors

        //- Construct from components
        NoMppicDamping
        (
            solidParticleCouplingCloud& spc,
            const dictionary& dict
            
        );

        //- Construct copy
//         NoDamping(const NoDamping<CloudType>& cm);

        //- Construct and return a clone
//         virtual autoPtr<DampingModel<CloudType> > clone() const
//         {
//             return autoPtr<DampingModel<CloudType> >
//             (
//                 new NoDamping<CloudType>(*this)
//             );
//         }


    //- Destructor
    virtual ~NoMppicDamping();


    //- Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            solidParticleCoupling& p,
            const scalar deltaT
        ) const;

        //- Return the model 'active' status
        virtual bool active() const;
        
        virtual void cacheFields(const bool store);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace DampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "NoDamping.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
