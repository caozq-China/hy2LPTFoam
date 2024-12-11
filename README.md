hy2LPTFoam is an open-source code for solving non-equilibrium high-speed gas-paricle flows.
For more detail, please find the paper: https://doi.org/10.3390/aerospace11090742

hy2LPTFoam is based on hy2Foam in hyStrath:https://github.com/hystrath/hyStrath.

Please install hy2Foam with OpenFOAM-v1706 before compile any code of hy2LPTFoam.

Once you have installed hy2Foam, you can copy the folder to your working space, e.g. XXXX-v1706.
The code can be compiled in the following steps:
1.Go to src/lagrangian/basic, run "wmake libso" in the terminal
2.Go to src/lagrangian/LPT, run "wmake libso" in the terminal
3.Go to applications/preProcessing/solidInitialise, copy BCs/ in applications/solvers/hy2LPTFoam to current directory and run "wmake libso" in BCs, then run "wmake" in the terminal to compile solidInitialise
4.Go to applications/solvers/hy2LPTFoam/BCs, run "wmake libso" in the terminal
5.Go to applications/solvers/hy2LPTFoam, run "wmake" in the terminal

A tutorial case corresponding to "MSRO body:hypersonic non-equilibrium flow during Mars entry" with particle mass fraction 0.014% is uploaded.
