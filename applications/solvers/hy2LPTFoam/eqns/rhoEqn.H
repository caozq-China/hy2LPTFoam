
// --- Solve density
//- phi is defined in createThermoFields.H ,equals to mesh.Sf() & fvc::interpolate(rhoU)
solve(fvm::ddt(thetaC,rho) + fvc::div(thetaCfPhi));

