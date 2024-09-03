
// --- Solve vibrational energy equation, INVISCID
thermalRelaxModelVT().correct();
volScalarField QVT = thermalRelaxModelVT().QVT();

if (containsElectrons)
{
    thermalRelaxModelHE().correct();
    QVT += thermalRelaxModelHE().QHE();
}

fvScalarMatrix rhoEvInviscid
(
    fvm::ddt(thetaC,rhoEv)
  + fvc::div(thetaCf*phiEvk[0])
//   - QVT*thetaC
  - QVT
  - reaction->Scv()//- chemistry-vibrational coupling source term
//     fvm::ddt(rhoEv)
//   + fvc::div(phiEvk[0])
//   - QVT*thetaC
//   - reaction->Scv()//- chemistry-vibrational coupling source term
);

if (containsElectrons)
{
    //- thermo.pe -- Electron pressure
    rhoEvInviscid += thermo.pe()*fvc::div(U) - reaction->Siir();
}

rhoEvInviscid.solve
(
    mesh.solver("rhoEv")
);

ev.correctBoundaryConditions();
rhoEv.boundaryFieldRef() = rho.boundaryField()*ev.boundaryField();

ev = rhoEv/rho;

