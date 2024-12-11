
p.ref() = rho() / psi();//- psi() = 1/(RT)
p.correctBoundaryConditions();
rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField();
