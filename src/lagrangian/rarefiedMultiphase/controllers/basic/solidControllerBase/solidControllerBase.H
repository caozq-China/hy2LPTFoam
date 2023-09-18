/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::solidControllerBase

Description
    Basic/abstract class of a controller

SourceFiles
    solidControllerBase.C

\*---------------------------------------------------------------------------*/

#ifndef solidControllerBase_H
#define solidControllerBase_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "vectorField.H"
#include "tensorField.H"
#include "Random.H"
#include "timeData.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParticleCouplingCloud;
class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class solidControllerBase Declaration
\*---------------------------------------------------------------------------*/

class solidControllerBase
{
protected:

    // Protected data

        const polyMesh& mesh_;

        solidParticleCouplingCloud& cloud_;

        Random& rndGenS_;

        //- Properties dictionary
        dictionary controllerDict_;

        dictionary timeDict_;

        timeData timeData_;

        //- name of face zone
        word regionName_;
        label regionId_;


        bool control_;
        bool readStateFromFile_;


        // set all the properties below from model if required

            //- target single valued controller
            bool singleValueController_;
            
            vector velocity_;
            scalar density_;
            scalar temperature_;
//             scalar pressure_;

//             tensor strainRate_;
//             vector tempGradient_;

            //- targeted field controller
            bool fieldController_;
            
            vectorField velocities_;
            scalarField densities_;
            scalarField temperatures_;
            scalarField pressures_;

            //- writing out properties
            bool writeInTimeDir_;
            bool writeInCase_;

            template<class Type>
            Type averageProperty
            (
                const Type& value,
                const Field<Type>& fieldValues
            ) const;


public:

    //- Construct from components
    solidControllerBase
    (
        const polyMesh& mesh,
        solidParticleCouplingCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~solidControllerBase() = default;


    // Member Functions

        void updateTime();

        void writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const List<Pair<scalar>>& data
        );

        //- scalar field
        void writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const scalarField& xAxis,
            const scalarField& yAxis
        );

        //- vector field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData
        );

        //- tensor field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const tensorField& yData
        );

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- calculate any required properties
        virtual void calculateProperties() = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);


        // Access

            //- Return the control elements
            virtual const labelList& controlZone() const = 0;

            //- return the control zone name
            virtual const word& regionName() const;

            //- return the targeted values
            virtual scalar density() const;
            virtual scalar& density();

            virtual const vector& velocity() const;
            virtual vector& velocity();

            virtual scalar temperature() const;
            virtual scalar& temperature();

//             virtual scalar pressure() const;
//             virtual scalar& pressure();
/*
            virtual const tensor& strainRate() const;
            virtual tensor& strainRate();*/
/*
            virtual const vector& tempGradient() const;
            virtual vector& tempGradient();*/

            //- return the targeted fields
            virtual const scalarField& densityField() const;
            virtual scalarField& densityField();

            virtual const vectorField& velocityField() const;
            virtual vectorField& velocityField();

            virtual const scalarField& temperatureField() const;
            virtual scalarField& temperatureField();

//             virtual const scalarField& pressureField() const;
//             virtual scalarField& pressureField();

            virtual bool singleValueController() const;
            virtual bool& singleValueController();

            virtual bool fieldController() const;
            virtual bool& fieldController();

            virtual bool writeInTimeDir() const;
            virtual bool writeInCase() const;


            virtual scalar avReqDensity();
            virtual vector avReqVelocity();
            virtual scalar avReqTemperature();
//             virtual scalar avReqPressure();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidControllerBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
