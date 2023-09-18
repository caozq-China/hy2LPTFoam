/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::dsmcFaceTracker

Description

SourceFiles
    dsmcFaceTrackerI.H
    dsmcFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcFaceTrackerFO_H
#define dsmcFaceTrackerFO_H

#include "vector.H"
#include "fvMesh.H"
#include "dsmcParcel.H"
#include "surfaceFields.H"
#include "CloudFunctionObject.H"
#include "CloudWithModels.H"
#include "dsmcCloud.H"

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dsmcFaceTracker Declaration
\*---------------------------------------------------------------------------*/

class dsmcFaceTrackerFO
:
    public CloudFunctionObject<dsmcCloud>
{
    // Private data

        List<scalarField> parcelIdFlux_;
        List<scalarField> massIdFlux_;


public:

    TypeName("dsmcFaceTrackerFO");

    // Constructors

        //- Construct from dictionary
        dsmcFaceTrackerFO
        (
            const dictionary& dict,
            dsmcCloud& owner,
            const word& modelName
        );

        //- Construct copy
        dsmcFaceTrackerFO(const dsmcFaceTrackerFO& ft);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<dsmcCloud>> clone() const
        {
            return autoPtr<CloudFunctionObject<dsmcCloud>>
            (
                new dsmcFaceTrackerFO(*this)
            );
        }

    //- Destructor
    virtual ~dsmcFaceTrackerFO() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve();

        //- Post-face hook
        virtual void postFace(const dsmcParcel& p, bool& keepParticle);

        inline const List<scalarField>& parcelIdFlux() const;

        inline const List<scalarField>& massIdFlux() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcFaceTrackerFOI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
