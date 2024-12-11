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

Description

\*---------------------------------------------------------------------------*/

#include "solidFieldProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidFieldProperties::solidFieldProperties
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    solidFieldPropertiesDict_
    (
        IOobject
        (
            "solidFieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(),
    fieldNames_(),
    fieldIds_(),
    fields_()
{}


solidFieldProperties::solidFieldProperties
(
    Time& t,
    const polyMesh& mesh,
    solidParcelCloud& cloud
)
:
    time_(t),
    solidFieldPropertiesDict_
    (
        IOobject
        (
            "solidFieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(solidFieldPropertiesDict_.lookup("solidFields")),
    fieldNames_(fieldList_.size()),
    fieldIds_(fieldList_.size()),
    fields_(fieldList_.size())
{
    if(fields_.size() > 0 )
    {
        Info << "Creating solid phase fields: " << nl << endl;
        
        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();

            fields_[f] = solidField::New(time_, mesh, cloud, fieldIDict);

            fieldNames_[f] = fields_[f]->type();
            fieldIds_[f] = f;
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void solidFieldProperties::createFields()
{
    Info << nl << "Initialising the measurement solid fields" << nl << endl;

    forAll(fields_, f)
    {
        fields_[f]->createField();
    }
}


//- this function is to be called at the beginning of the MD time-step,
//  in order to update the time scheme used by each measurement property.
void solidFieldProperties::updateTimeInfo()
{

    forAll(fields_, f)
    {
        fields_[f]->updateTime();
    }
}


void solidFieldProperties::calculateFields()
{

    forAll(fields_, f)
    {
        fields_[f]->calculateField();
    }
}



//- Note, not all fields automatically write out to hard disc. 
void solidFieldProperties::writeFields()
{
    const Time& runTime = time_;

    fileName timePath(runTime.path()/runTime.timeName()/"uniform");

    if (runTime.writeTime())
    {
        if (Pstream::master())
        {
            if (!isDir(timePath))
            {
                mkDir(timePath);
            }
        }
    }

    forAll(fields_, f)
    {
        fields_[f]->timePath() = timePath;
        fields_[f]->writeField();
    }
   
    if(runTime.writeTime())
    {
        //- Checking for modifications in the IOdictionary
        //  this allows for run-time tuning of any parameters.

        // NOTES: 
        // At the moment, the dictionary is forced to be re-read every write-interval, 
        // and properties within the abstract and models are re-set to. 
        // The "ideal" case is to have the code identify when the dictionary has been
        // modified, before re-reading it in again. Unfortunately the .modified() function
        //  is not working properly. 

        fieldList_.clear();
    
        fieldList_ = solidFieldPropertiesDict_.lookup("solidFields");
    
        forAll(fields_, f)
        {
            const entry& fieldEntry = fieldList_[f];
            fields_[f]->updateProperties(fieldEntry.dict());
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
