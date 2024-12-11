/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

Application
    dsmcFoamInitialise

Description
    Initialise a case for dsmcFoam by reading the initialisation dictionary
    system/dsmcInitialise

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solidParcelCloud.H"
#include "rho2HTCModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOdictionary solidInitialiseDict
    (
        IOobject
        (
            "solidInitialiseDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Creating reaction model\n" << endl;

    autoPtr<hTC2Models::rho2HTCModel> reaction
    (
        hTC2Models::rho2HTCModel::New(mesh)
    );

    rho2ReactionThermo& thermo = reaction->thermo();

    thermo.validate(args.executable(), "h", "e");

    //---------------------------------------------------------
    // Initialisation of the two-temperature model
    //---------------------------------------------------------

    thermo.initialise2T();

    // volScalarField rhoG
    // (
    //     IOobject
    //     (
    //         "rho",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    // );

    volScalarField Mach
    (
        IOobject
        (
            "Mach",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Mach", dimless, 0)
    );

    volScalarField specificHeatRatio
    (
        IOobject
        (
            "specificHeatRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("specificHeatRatio", dimless, 0)
    );

    volScalarField viscosityIndex
    (
        IOobject
        (
            "viscosityIndex",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("viscosityIndex", dimless, 0)
    );

    volVectorField UG
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("U",dimensionSet(0,1,-1,0,0,0,0), Zero)
    );

    // --- Initialisation of lagrangian fields
    #include "readGravitationalAcceleration.H"

    Info<< "Initialising solid particles for Time = " << runTime.timeName() << nl <<endl;

    solidParcelCloud spc
    (
        runTime,
        "spc",
        mesh,
        // rhoG,
        UG,
        Mach,
        viscosityIndex,
        specificHeatRatio,
        g,
        thermo,
        solidInitialiseDict,
        true
    );

//     label totalMolecules = dsmc.size();
//
//     if (Pstream::parRun())
//     {
//         reduce(totalMolecules, sumOp<label>());
//     }
//
//     Info<< nl << "Total number of molecules added: " << totalMolecules
//         << nl << endl;

    spc.updateCellOccupancy();
    spc.updateTheta();
    
    IOstream::defaultPrecision(15);

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing solidParcelCloud."
            << nl << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
