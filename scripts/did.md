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

## 2020-08-17 23:32:28+02:00

 * test script `test_pce_mc.py`
 * pce very exact -- mc horribly slow convergence
 * but works well -- will include this in the pub

## 2020-08-25 10:45:31+02:00

 * TODO: rerun big sims to also have the variance
 * TODO: already prepare for code pub
 * [ ] `spacetime_galerkin_pod` to `pip`
 * [ ] add the pips for the dependencies
 * [ ] prop this to mechthild

# Adding some RB Functionality for the revision

## 2021-06-09 14:22:24+02:00

 * can use the MC routines to compute the truth snapshots 
 * TODO: check if the list represents the param list
 * DONE: it does!!!
 * basic RB functionality -- enriching the basis -- like working
   * TODO: need check with the errors and whether the right vecs are returned
   * DONE!

## 2021-07-26 08:10:18+02:00

 * RB works like very well now
 * NEXT: use the RB to comp Ey and so on

## 2022-03-01 21:03:50+01:00

 * Will use the RB with PCE like with genpod
 * DONE: compute an M-ortho basis of every RB 
    * to do the projection
    * to be inline with genpod
 * TODO: clean up code

## 2022-03-03 20:51:07+01:00

 * DONE: do some caching
   * save `expvals` and `variances` for all `(N, nua, nub)`
   * define `pcesnapshotdim` separately (independent from `pceplease`)
   * one big json...
     * `N`
       * `nua-nub`
         * `uniform`
           * `pcedim`
             * `expv`
             * `vrnc`

## 2022-03-23 14:52:35+01:00

 * TODO: define/comp the `rb` experiments
 * TODO: add support for `beta-dist` see Notion note for the formulas
