
// --- Solve momentum equation, INVISCID
solve(fvm::ddt(thetaC,rhoU) + fvc::div(thetaCf*phiUp));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
if (hyDev)
{
    thermo.limitVelocityAtWallBoundary(U.boundaryFieldRef());
}

rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();


// --- Solve momentum equation, VISCOUS
if (!inviscid)
{

    fvVectorMatrix UEqn
    (
        fvm::ddt(thetaC, rho, U) - fvc::ddt(thetaC, rho, U)
      - fvm::laplacian(thetaC*muEff, U)
      - fvc::div(thetaC*tauMC)
     ==
      mhd->lorentzForce()
      + fvOptions(rho, U)
      + particles.SU(U)

    );

    UEqn.relax();
    
    fvOptions.constrain(UEqn);
    
    UEqn.solve();
    
    fvOptions.correct(U);
    
    rhoU = rho*U;
}
