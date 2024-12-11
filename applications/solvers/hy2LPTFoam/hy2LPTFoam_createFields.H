#include "createTimeControls.H"
#include "createThermoFields.H"
#include "readTimeControls.H"

// --- Initialisation of reacting and two-temperature fields
#include "createReactingFields.H"
#include "createVibrationalFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info << "Creating volume fraction field \n" << endl;



// --- Initialisation of lagrangian fields
#include "readGravitationalAcceleration.H"

bool evolveCloud = Switch(runTime.controlDict().lookup("evolveParticles"));

Info<< "Constructing Lagrangain Cloud " << endl;
Info<<"g="<<g<<endl;
solidParcelCloud particles
(
    runTime,
    "spc",
    mesh,
    U,
    Mach,
    specificHeatRatio,
    viscosityIndex,
    g,
    thermo
);


//- create volume fraction
// Particle fraction upper limit
scalar thetaCMin
(
    1.0
- readScalar
    (
        particles.particleProperties().subDict("constantProperties")
    .lookup("thetaMax")
    )
);
Info<<"thetaCMin="<<thetaCMin<<endl;
volScalarField thetaC
(
    IOobject
    (
        "thetaC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

particles.updateTheta();
thetaC = max(1.0 - particles.theta(), thetaCMin);
thetaC.correctBoundaryConditions();

volScalarField thetaCfilter
(
    IOobject
    (
        "thetaCfilter",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thetaC
);

surfaceScalarField thetaCf("thetaCf", fvc::interpolate(thetaCfilter));
//- phi is defined as "phi("phi", mesh.Sf() & fvc::interpolate(rhoU));"
surfaceScalarField thetaCfPhi(thetaCf*phi);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// --- Initialisation of MHD fields
#include "createEMFields.H"
#include "createFvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "numerics/readFluxScheme.H"
// --- Write fields in the 0 folder
#include "write/write0.H"

dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

// --- Upwind interpolation of primitive fields on faces
#include "numerics/upwindInterpolation.H"
#include "numerics/fluxesCalculation.H"

// --- Time control and local time stepping
#include "numerics/compressibleCourantNo.H"
#include "LTS/setInitialrDeltaT.H"
#include "setInitialDeltaT.H"
