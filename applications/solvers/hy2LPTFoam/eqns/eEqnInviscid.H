
// --- Solve energy, INVISCID
surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
      + fvc::dotInterpolate(mesh.Sf(), tauMC)
    )
  & aByU
);

solve
(
    fvm::ddt(thetaC,rhoE)
    + fvc::div(thetaCf*phiEp)//- correspending to rhoEFlux
    - fvc::div(thetaCf*sigmaDotU)
  ==
      reaction->Sh()
    + mhd->jouleHeating(U)
);

e = rhoE/rho - 0.5*magSqr(U);

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

rhoE.boundaryFieldRef() =
    rho.boundaryField()
  * (
        e.boundaryField() + 0.5*magSqr(U.boundaryField())
    );

