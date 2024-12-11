/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    solidVolFields

Description
    Geometric field for temperature, density, velocity, pressure.



SourceFiles
    solidVolFieldsI.H
    solidVolFields.C
    solidVolFieldsIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidVolFields_H
#define solidVolFields_H

#include "solidField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "solidParcelCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidVolFields Declaration
\*---------------------------------------------------------------------------*/

class solidVolFields
:
    public solidField
{

private:

    // Private data
        dictionary propsDict_;
       
        //- Wall unit vectors - temporary storage, should be done in 
        // boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

        label sampleInterval_;
        label sampleCounter_;
        scalar nTimeSteps_;
        word fieldName_;

        labelList typeIds_;

        volScalarField nParcelMean_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField q_;
        volVectorField UMean_;

        scalarField nParcelsCum_;
        scalarField nCum_;
        scalarField mCum_;
        vectorField momentumCum_;

        List<labelList> boundaryCells_;
        
        // boundary measurements
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> qBF_;
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;
        
        
        bool averagingAcrossManyRuns_;
        
    // Private Member Functions

        void readIn();

        void writeOut();
       
        void calculateWallUnitVectors();


public:

    //- Runtime type information
    TypeName("solidVolFields");

    // Constructors

        //- Construct from components
    solidVolFields
    (
        Time& t,
        const polyMesh& mesh,
        solidParcelCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    virtual ~solidVolFields() = default;


    // Member Functions

    //- create initial field
    virtual void createField();

    //- calculate field
    virtual void calculateField();

    //- write field
    virtual void writeField();
    
    volScalarField& surfaceHeatFlux();
    
    const word& fieldName() const;

    virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

