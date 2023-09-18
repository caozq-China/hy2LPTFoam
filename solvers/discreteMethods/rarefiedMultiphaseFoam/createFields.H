Info<< nl << "Constructing dsmcCloud " << endl;

dsmcCloud dsmc(runTime, "dsmc", mesh); 

Info<< nl << "Constructing solidParticleCouplingCloud " << endl;

solidParticleCouplingCloud spc(runTime, "spc", mesh,  &dsmc);
