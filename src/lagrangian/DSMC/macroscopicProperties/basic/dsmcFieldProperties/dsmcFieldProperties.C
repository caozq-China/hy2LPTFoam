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

\*---------------------------------------------------------------------------*/

#include "dsmcFieldProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcFieldProperties::dsmcFieldProperties
(
    const Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    dsmcFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(),
    fieldNames_(),
    fieldIds_(),
    fields_(),
    dsmcWeight_()
{}


Foam::dsmcFieldProperties::dsmcFieldProperties
(
    const Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    dsmcFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(dsmcFieldPropertiesDict_.lookup("dsmcFields")),
    fieldNames_(fieldList_.size()),
    fieldIds_(fieldList_.size()),
    fields_(fieldList_.size()),
    dsmcWeight_(mesh.nCells(), 0.0)
{
    if (fields_.size() > 0)
    {
        Info << "Creating fields: " << nl << endl;

        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();

            fields_[f] = dsmcField::New(time_, mesh, cloud, fieldIDict);

            fieldNames_[f] = fields_[f]->type();
            fieldIds_[f] = f;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::dsmcFieldProperties::createFields()
{
    Info << nl << "Initialising the measurement fields" << nl << endl;

    forAll(fields_, f)
    {
        fields_[f]->createField();
    }
}


void Foam::dsmcFieldProperties::updateTimeInfo()
{
    // This function is to be called at the beginning of the DSMC time-step,
    // in order to update the time scheme used by each measurement property.
    forAll(fields_, f)
    {
        fields_[f]->updateTime();
    }
}


void Foam::dsmcFieldProperties::calculateFields()
{
    forAll(fields_, f)
    {
        fields_[f]->calculateField();
    }
}


void Foam::dsmcFieldProperties::writeFields()
{
    // Note - not all fields automatically write out to file
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

    forAll(fields_, i)
    {
        fields_[i]->timePath() = timePath;
        fields_[i]->writeField();
    }

    if (runTime.writeTime())
    {
        // Checking for modifications in the IOdictionary
        // this allows for run - time tuning of any parameters.

        // NOTES:
        // At the moment, the dictionary is forced to be re-read every
        // write-interval, and properties within the abstract and models are
        // re-set to zero.
        // The "ideal" case is to have the code identify when the dictionary has
        // been modified, before re-reading it in again. Unfortunately the
        // .modified() function is not working properly.

        fieldList_.clear();

        fieldList_ = dsmcFieldPropertiesDict_.lookup("dsmcFields");

        forAll(fields_, i)
        {
            const entry& fieldEntry = fieldList_[i];
            fields_[i]->updateProperties(fieldEntry.dict());
        }
    }
}




// ************************************************************************* //
