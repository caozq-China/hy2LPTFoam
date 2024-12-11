// --- Re-read run-time editable dictionaries
#include "runTimeEditing/hTCPropertiesDictModification.H"
#include "runTimeEditing/transportDictModification.H"
#include "runTimeEditing/twoTemperatureDictModification.H"

// --- Upwind interpolation of primitive fields on faces
#include "numerics/upwindInterpolation.H"
// --- Compute numerical fluxes
#include "numerics/fluxesCalculation.H"

// --- Time control and local time stepping
#include "numerics/compressibleCourantNo.H"
#include "readTimeControls.H"
#include "LTS/readLTSTimeControls.H"
#include "setDeltaT.H"

runTime++;
Info<<nl
    <<"========================================"<<nl
    << "Time = " << runTime.timeName() << nl << endl;

volScalarField muEff(turbulence->muEff());
volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

/*if (hyDev)
{
    Zrot = ZrotInf/(1.0 + sqr(pi)/2.0*sqrt(TZ/thermo.T())
        + (sqr(pi)/4.0 + pi)*TZ/thermo.T());
    Z = 3.0*pi/(4.0*(thermo.zetar()+3.0))*Zrot;

    mub = mu*(2.0*thermo.zetar()*Z)/(3.0*(thermo.zetar()+3.0));
    bulkToShearRatio = mub/mu;

    tauMC = tauMC + mub*tr(fvc::grad(U))*I;
}*/

/*---------------------START EVOLVE PARTICLES------------------------*/
if(evolveCloud)
{
    //- update gas volume fraction, cell viscosity index, Mach number 
    //- update cell viscosity index
    forAll(rho.primitiveField(),celli)
    {
        //- reset
        label nSpecies = 0;
        viscosityIndex[celli] = 0;

        forAll(Y,specie)
        {
            //- acquire viscosity index of a specie
            const scalar viscoIndex = composition.omega(specie);
            if(Y[specie][celli]>VSMALL)
            {
                viscosityIndex[celli] += viscoIndex;
                nSpecies++;
            }
            
        }
        if(nSpecies>0)
        {
            viscosityIndex[celli] /= nSpecies;
        }
    }

    //-already tested viscosityIndex works.
    Mach = mag(U)/sqrt(thermo.Cp_t()/thermo.Cv_t()/psi);
    specificHeatRatio = thermo.Cp_t()/thermo.Cv_t();
    // End of update gas properties

    if(particles.solution().useGasFilter())
    {
        particles.filteringGasProperties();
    }
    else
    {
        particles.updateGasProperties();
    }

    particles.evolve();//- include updateTheta

    if(particles.solution().smoothTheta())
    {
        particles.DBS().diffusion(particles.theta());
    }
    
    //- update continuous phase volume fraction fields
    thetaC = max(1-particles.theta(), thetaCMin);
    thetaC.correctBoundaryConditions();

    thetaCfilter = thetaC;
    //- acquire the continuous phase volume fraction at cell surfaces
    thetaCf = fvc::interpolate(thetaCfilter);
    thetaCfPhi = thetaCf*phi;  

    if(particles.solution().smoothSourceTerms())
    {
        particles.smoothSourceTerms();
    }  
}

// particles.smoothSourceTerms();
// fvVectorMatrix cloudSU(particles.SU());
/*---------------------END EVOLVE PARTICLES------------------------*/

// --- Local time stepping (LTS)
if (activateLTS)
{
    #include "LTS/setrDeltaT.H"
}

// --- Re-set the switch/counter that serves as a warning if the
//     temperature goes unbounded
thermo.resetTemperatureBoundingInfo();

Info<<"- Solving CFD equations -"<<endl;
// --- Solve continuity equation
#include "eqns/rhoEqn.H"

// --- Solve species transport and reaction equations
#include "eqns/YEqn.H"

// --- Solve momentum equations
#include HY2FOAM_EXTERNAL_FILE_MOMENTUM_EQNS

if (downgradeSingleT)
{
    // --- Solve the total energy equation
    //- inviscid
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else if (downgradeSingleTv)
{
    // --- Solve the vibrational energy equation
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/evEqnInviscid_SingleTv.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include "eqns/evEqnViscous_SingleTv.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else if (downgradeSingleVibMode)
{
    // --- Solve the vibrational energy equations
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/evEqnInviscid.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include "eqns/evEqnViscous.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else
{
    // ABORTIVE WORK
    // --- Solve the vibrational energy equations (one per vib. mode)
    // --- Solve the total energy equation
    //- inviscid
    //#include "eqns/evEqn_MultiModes.H"
    //#include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    //#include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}

// --- Electro-magnetic fields calculation
#include "eqns/BEqn.H"

// --- Pressure field calculation
#include "eqns/pEqn.H"

// --- Enable input from / coupling with another application
#ifdef HY2FOAM_EXTERNAL_FILE_HYBRID_COUPLING
#include HY2FOAM_EXTERNAL_FILE_HYBRID_COUPLING
#endif

turbulence->correct();

mhd->update(U);

rarefactionParameters().correct(U);


// --- Print a report in the log file if temperature had to be bounded
thermo.temperatureBoundingReport();

if (runTime.outputTime())//- output of outputTime() is writeTime_ 
{
    runTime.write();
    #include "write/write.H"

    #ifdef HY2FOAM_EXTERNAL_FILE_OUTPUT
    #include HY2FOAM_EXTERNAL_FILE_OUTPUT
    #endif
}

previousIterationTime =
    max(runTime.elapsedCpuTime()-currentIterationTime, 1e-3);

Info<< "Phase no " << noRestart << "." << noSubRestart
    << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << "  Iteration no " << noIteration<<" (" << previousIterationTime
    << " s)" << nl << endl;

currentIterationTime = runTime.elapsedCpuTime();
noIteration += 1;
