import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

import spacetime_galerkin_pod.chaos_expansion_utils as ceu
import spacetime_galerkin_pod.ten_sor_utils as tsu
import spacetime_galerkin_pod.gen_pod_utils as gpu
from spacetime_galerkin_pod.ldfnp_ext_cholmod import SparseFactorMassmat

import gen_pod_uq.mc_pce_utils as mpu

from circle_subsec import get_problem


get_linop, get_sol, get_output, problemfems = get_problem()
print(problemfems['mmat'].shape[0])

uncdims = 5

basenu = 1e-3
varia = 0.
varib = 5e-4
nua, nub = basenu+varia, basenu+varib

mcits, mcruns = 6, 100  # 200
# pcedimlist = [2, 3, 5]
pcedimlist = [4]  # , 3, 4, 5]  # , 7]

mcplease = False
pceplease = False
plotplease = False
# ## make it come true
# mcplease = True
pceplease = True
# plotplease = True

basenulist = [basenu]*uncdims
basey = get_output(basenulist)
print('y(estxnu)={0}'.format(basey))
# import ipdb
# ipdb.set_trace()


# ## CHAP Monte Carlo
if mcplease:
    varinu = basenu + (varib-varia)*np.random.rand(mcits*mcruns, uncdims)
    expvnu = np.average(varinu, axis=0)
    print('expected value of nu: ', expvnu)
    varinulist = varinu.tolist()
    mcout, expvnu = mpu.run_mc_sim(varinulist, get_output, verbose=True)

    curyplotfignum = 101 if plotplease else None
    cury = get_output(expvnu.tolist(), plotfignum=curyplotfignum)
    print('y(estxnu)={0}'.format(cury))

    if plotplease:
        plt.figure(89)
        plt.plot(mcout, '.')
        plt.show()


# ## CHAP Polynomial Chaos Expansion
if pceplease:
    for pcedim in pcedimlist:
        abscissae, weights, compexpv = mpu.\
            setup_pce(distribution='uniform',
                      distrpars=dict(a=nua, b=nub),
                      pcedim=pcedim, uncdims=uncdims)
        # abscarray, weightsarray = np.array(abscissae), np.array(weights)
        ysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
                                             uncdims=uncdims,
                                             abscissae=abscissae)
        expy = compexpv(ysoltens)
        print('PCE({0}): E(y): {1}'.format(pcedim, expy))

# ## CHAP genpod
pcedim = pcedimlist[0]
mmat = problemfems['mmat']

abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)
facmy = SparseFactorMassmat(mmat)
pcemmat = sps.csc_matrix(sps.dia_matrix((weights, 0), shape=(pcedim, pcedim)))
facmpce = SparseFactorMassmat(pcemmat)

basisfrom = 'pce'
basisfrom = 'mc'
poddimlist = [5, 10, 20]  # , 40]
nmcsnapshots = 5*pcedim**uncdims

pcewmat = sps.dia_matrix((weights, 0), shape=(pcedim, pcedim))
pcewmatfac = sps.dia_matrix((np.sqrt(weights), 0), shape=(pcedim, pcedim))

mfl = [facmy.F]
mfl.extend([pcewmatfac]*uncdims)

if basisfrom == 'pce':
    ysoltens = mpu.run_pce_sim_separable(solfunc=get_sol,
                                         uncdims=uncdims,
                                         abscissae=abscissae)

    def get_pod_vecs(poddim=None):
        return tsu.modeone_massmats_svd(ysoltens, mfl, poddim)


elif basisfrom == 'mc':
    varinu = basenu + (varib-varia)*np.random.rand(nmcsnapshots, uncdims)
    expvnu = np.average(varinu, axis=0)
    print('expected value of nu: ', expvnu)
    varinulist = varinu.tolist()
    mcout, expvnu = mpu.run_mc_sim(varinulist, get_sol, verbose=True)
    pceymat = np.array(mcout).T
    lypceymat = facmy.Ft*pceymat

    def get_pod_vecs(poddim=None):
        ypodvecs = gpu.get_ksvvecs(sol=lypceymat, poddim=poddim,
                                   plotsvs=plotplease, labl='Singular Values')
        return ypodvecs

nulist = [basenu]*5
nuarray = np.array(nulist)

# lypceymat = pceymat
for poddim in poddimlist:

    ypodvecs = get_pod_vecs(poddim)
    lyitVy = facmy.solve_Ft(ypodvecs)

    def red_out_func(parlist):
        return get_output(parlist, podmat=lyitVy)

    if plotplease:
        yfull = get_output(nuarray.tolist(), plotfignum=222)
        yred = get_output(nuarray.tolist(), plotfignum=111, podmat=lyitVy)

    if pceplease:
        for pcedim in pcedimlist:
            abscissae, weights, compredexpv = mpu.\
                setup_pce(distribution='uniform',
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcedim, uncdims=uncdims)
            redysoltens = mpu.run_pce_sim_separable(solfunc=red_out_func,
                                                    uncdims=uncdims,
                                                    abscissae=abscissae)
            redexpy = compredexpv(redysoltens)
            print('pcedim={0:2.0f}, poddim={2:2.0f}, exypce={1}'.
                  format(pcedim, redexpy-expy,
                         poddim))

plt.show()
