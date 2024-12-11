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
    solidFieldProperties

Description

	Stores all the information for the solidFieldProperties models
    defined within the fieldPropertiesDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef solidFieldProperties_H
#define solidFieldProperties_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "solidField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidFieldProperties Declaration
\*---------------------------------------------------------------------------*/

class solidFieldProperties
{

    // Private data

        //- The entire field-properties dictionary 
        //  (it contains multiple subDictionaries)

        Time& time_;

        IOdictionary solidFieldPropertiesDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List< autoPtr<solidField> > fields_;

public:

    // Constructors

        //- Null Constructor (mdInitialise)
        solidFieldProperties
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFoam
        solidFieldProperties
        (
            Time& t,
            const polyMesh& mesh,
            solidParcelCloud& cloud
        );

        virtual ~solidFieldProperties() = default;

    // Member Functions


        void updateTimeInfo();

        void createFields();

        void calculateFields();

        void writeFields();

        inline List< autoPtr<solidField> >& fields();
        inline const List< autoPtr<solidField> >& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;

        inline volScalarField& surfaceHeatFlux();

    // Ostream Operator
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidFieldPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
