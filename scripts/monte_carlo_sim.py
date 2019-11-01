import itertools

import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

import spacetime_galerkin_pod.chaos_expansion_utils as ceu
import spacetime_galerkin_pod.ten_sor_utils as tsu
import spacetime_galerkin_pod.gen_pod_utils as gpu

from spacetime_galerkin_pod.ldfnp_ext_cholmod import SparseFactorMassmat

from circle_subsec import get_problem

get_linop, get_sol, get_output, problemfems = get_problem()
print(problemfems['mmat'].shape[0])

uncdims = 5

basenu = 1e-3
varia = 0.
varib = 5e-4
nua, nub = basenu+varia, basenu+varib

mcits, mcruns = 20, 5000  # 200
# pcedimlist = [2, 3, 5]
pcedimlist = [5]  # , 7]

mcplease = False
pceplease = False
plotplease = False
# ## make it come true
# mcplease = True
pceplease = True
plotplease = True

# ## CHAP Monte Carlo
if mcplease:
    varinu = basenu + (varib-varia)*np.random.rand(mcits*mcruns, uncdims)
    expvnu = np.average(varinu, axis=0)
    print('expected value of nu: ', expvnu)
    varinulst = []
    for uncdim in range(uncdims):
        varinulst.append(varinu[:, uncdim].tolist())

    # estxnu, estxy = 0, 0
    nulist = [basenu]*5
    ylist = []

    for mitk in range(mcits):
        for mck in range(mcruns):
            for uncdim in range(uncdims):
                nulist[uncdim] = varinulst[uncdim].pop(0)
            cury = get_output(nulist, plotfignum=None)
            ylist.append(cury)
        estxy = np.average(np.array(ylist))
        print('mc:{0}/{1}: estxy={2}'.
              format((mitk+1)*mcruns, mcits*mcruns, estxy))

    for uncdim in range(uncdims):
        nulist[uncdim] = expvnu[uncdim]

    curyplotfignum = 101 if plotplease else None
    cury = get_output(nulist, plotfignum=curyplotfignum)
    print('y(estxnu)={0}'.format(cury))

    if plotplease:
        plt.figure(89)
        plt.plot(ylist, '.')
        # plt.figure(98)
        # plt.plot(xnulist, '.')
        plt.show()


# ## CHAP Polynomial Chaos Expansion

if pceplease:
    ydim = 1  # dimension of the output
    for pcedim in pcedimlist:
        pceylist = []
        ypcedims = [ydim]
        ypcedims.extend([pcedim]*uncdims)
        ypcedims = tuple(ypcedims)

        abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)
        # abscarray, weightsarray = np.array(abscissae), np.array(weights)

        nulist = [basenu]*5
        nuarray = np.array(nulist)
        for idxtuple in itertools.product(np.arange(pcedim), repeat=uncdims):
            idxarray = np.array(idxtuple)
            nuarray[:uncdims] = abscissae[idxarray]
            pceylist.append(get_output(nuarray.tolist(), plotfignum=None))

        yrslttns = np.array(pceylist).reshape(ypcedims)
        yrslttns = tsu.tnsrtrnsps(yrslttns)
        exypce = 0
        # for kk, cw in enumerate(weights):
        for idxtuple in itertools.product(np.arange(pcedim), repeat=uncdims):
            idxarray = np.array(idxtuple)
            cw = (weights[idxarray]).prod()
            exypce += cw*yrslttns[idxtuple]

        print('pcedim={0:2.0f}, exypce={1}'.
              format(pcedim, 1./((varib-varia)**uncdims)*exypce))
        pceexy = 1./((varib-varia)**uncdims)*exypce

# ## CHAP genpod
pcedim = pcedimlist[0]
mmat = problemfems['mmat']
ydim = mmat.shape[0]
ypcedims = [ydim]
ypcedims.extend([pcedim]*uncdims)
ypcedims = tuple(ypcedims)
nulist = [basenu]*5
nuarray = np.array(nulist)
abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)

pceylist = []
for idxtuple in itertools.product(np.arange(pcedim), repeat=uncdims):
    idxarray = np.array(idxtuple)
    nuarray[:uncdims] = abscissae[idxarray]
    pceylist.append(get_sol(nuarray.tolist()))
yrslttns = np.array(pceylist).reshape(ypcedims)

facmy = SparseFactorMassmat(mmat)
pcemmat = sps.csc_matrix(sps.dia_matrix((weights, 0), shape=(pcedim, pcedim)))
facmpce = SparseFactorMassmat(pcemmat)

basisfrom = 'mc'
basisfrom = 'pce'
poddimlist = [5, 10, 20, 40]
nmcsnapshots = 5*pcedim**uncdims

if basisfrom == 'pce':
    pceymat = np.array(pceylist).T
    pceymat = pceymat[0, :, :]
elif basisfrom == 'mc':
    varinu = basenu + (varib-varia)*np.random.rand(nmcsnapshots, uncdims)
    # expvnu = np.average(varinu, axis=0)
    # print('expected value of nu: ', expvnu)
    varinulst = []
    for uncdim in range(uncdims):
        varinulst.append(varinu[:, uncdim].tolist())
    ylist = []
    for mck in range(nmcsnapshots):
        for uncdim in range(uncdims):
            nulist[uncdim] = varinulst[uncdim].pop(0)
        cury = get_sol(nulist)
        ylist.append(cury)
    pceymat = np.array(ylist).T
    pceymat = pceymat[0, :, :]

lypceymat = facmy.Ft*pceymat
# lypceymat = pceymat
for poddim in poddimlist:
    ypodvecs = gpu.get_ksvvecs(sol=lypceymat, poddim=poddim,
                               plotsvs=plotplease, labl='Singular Values')

    # massfaclist = [facmy.F]
    # massfaclist.extend([facmpce.F]*uncdims)
    # ypodvecs = tsu.modeone_massmats_svd(yrslttns, massfaclist, kdim=poddim)

    lyitVy = facmy.solve_Ft(ypodvecs)
    # lyitVy = ypodvecs

    if plotplease:
        yfull = get_output(nuarray.tolist(), plotfignum=222)
        yred = get_output(nuarray.tolist(), plotfignum=111, podmat=lyitVy)

    if pceplease:
        ydim = 1  # dimension of the output
        for pcedim in pcedimlist:
            pceylist = []
            ypcedims = [ydim]
            ypcedims.extend([pcedim]*uncdims)
            ypcedims = tuple(ypcedims)

            abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim,
                                                         a=nua, b=nub)
            # abscarray, weightsarray = np.array(abscissae), np.array(weights)

            nulist = [basenu]*5
            nuarray = np.array(nulist)
            for idxtuple in itertools.product(np.arange(pcedim),
                                              repeat=uncdims):
                idxarray = np.array(idxtuple)
                nuarray[:uncdims] = abscissae[idxarray]
                pceylist.append(get_output(nuarray.tolist(), plotfignum=None,
                                           podmat=lyitVy))

            yrslttns = np.array(pceylist).reshape(ypcedims)
            yrslttns = tsu.tnsrtrnsps(yrslttns)
            exypce = 0
            # for kk, cw in enumerate(weights):
            for idxtuple in itertools.product(np.arange(pcedim),
                                              repeat=uncdims):
                idxarray = np.array(idxtuple)
                cw = (weights[idxarray]).prod()
                exypce += cw*yrslttns[idxtuple]

            print('pcedim={0:2.0f}, poddim={2:2.0f}, exypce={1}'.
                  format(pcedim, 1./((varib-varia)**uncdims)*exypce-pceexy,
                         poddim))

plt.show()
