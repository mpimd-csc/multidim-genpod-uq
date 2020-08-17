## 19-12-12

 * 3D works
 * mesh 3 -- POD bas for PCE(3)
 * very good results for redPCE(3)
 * similarly good for redPCE(5)
 * poddim 10 is worse than 5 but much worse than 20

It's good that although trained for PCE(3), podPCE(5) works as the abscissae are
very different.

 * TODO: check the error PCE(3)-PCE(5) and see whether pod helps here.
 * TODO: benchmark with MC sims
 * DONE: MC sims confirm the PCE expv

## 20-01-05

 * added convection in z-direc
 * less diffusion: `5e-4`, more variation: `+-1e-4`
 * new mesh (refined at the corners)
 * goal: better convergence both in PCE and podvecs
 * reasonable mesh: 6, gonna check against 7...

## 20-01-08

 * rewrote mesh to better control refinement
 * convergence of the results in the range of 1e-5 only for ML>10 = #dof~200.000
