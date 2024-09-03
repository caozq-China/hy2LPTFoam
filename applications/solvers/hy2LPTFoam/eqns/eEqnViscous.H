
// --- Solve energy, VISCOUS
if (!inviscid)
{
    fvScalarMatrix eEqnViscous
    (
        fvm::ddt(thetaC,rho,e) - fvc::ddt(thetaC, rho, e)
      - fvm::laplacian(thermo.kappa()/thermo.CvMix(), e)//- heat transfer
      ==
        fvOptions(rho, e)
      + particles.SE(e)
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        eEqnViscous += 
            fvc::div
            (
                thetaC*
                speciesDiffusion().multiSpeciesHeatSource(),
                "div(diffusiveFluxes)"
            );
    }
    
    eEqnViscous.relax();

    fvOptions.constrain(eEqnViscous);
    
    eEqnViscous.solve();

    fvOptions.correct(e);

    et.ref() = e() - ev();

    if (downgradeSingleT)
    {
        e.correctBoundaryConditions();
        et.boundaryFieldRef() = e.boundaryField() - ev.boundaryField();
    }
    else
    {
        et.correctBoundaryConditions();
        e.boundaryFieldRef() = et.boundaryField() + ev.boundaryField();
    }

    rhoE = rho*(e + 0.5*magSqr(U));
}

if (hyLight)
{
    thermo.correct2T_Light();
}
else
{
    thermo.correct2T();
    thermo.correctOverallTemperature();
}
