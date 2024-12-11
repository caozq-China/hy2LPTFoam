
// --- Solve vibrational energy equation, VISCOUS
if (!inviscid)
{
    fvScalarMatrix evEqnViscous
    (
        fvm::ddt(thetaC, rho, ev) - fvc::ddt(thetaC, rho, ev)
      - fvm::laplacian(transport().kappave()/thermo.CvvelMix(), ev)
      ==
        fvOptions(rho, ev)
    //     fvm::ddt(rho, ev) - fvc::ddt(rho, ev)
    //   - fvm::laplacian(transport().kappave()/thermo.CvvelMix(), ev)
    //   ==
    //     fvOptions(rho, ev)
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        evEqnViscous +=
            fvc::div
            (
                thetaC*
                speciesDiffusion().multiSpeciesVEHeatSource(),
                "div(diffusiveFluxes)"
            );
    }

    evEqnViscous.solve();

    ev.correctBoundaryConditions();
    rhoEv = rho*ev;
}

