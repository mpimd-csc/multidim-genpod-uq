import itertools

import numpy as np
import matplotlib.pyplot as plt

import spacetime_galerkin_pod.chaos_expansion_utils as ceu

from circle_subsec import get_problem

get_linop, get_sol, get_output, problemfems = get_problem()

uncdims = 5

basenu = 1e-3
varia = 0.
varib = 5e-4

mcits, mcruns = 5, 500  # 200
ydim = 1  # dimension of the output
pcedimlist = [2, 3]

plotplease = True
plotplease = False

# ## CHAP Monte Carlo

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
    print('mc:{0}/{1}: estxy={2}'.format((mitk+1)*mcruns, mcits*mcruns, estxy))

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

paxlist = np.arange(uncdims+1).tolist()
paxlist.append(paxlist.pop(0))


def tnsrtrnsps(X, times=1):
    '''transpose the tensor (by cycling the dimensions)'''
    for k in range(times):
        X = np.transpose(X, paxlist)
    return X


nua, nub = basenu+varia, basenu+varib

for pcedim in pcedimlist:
    pceylist = []
    ypcedims = [ydim]
    ypcedims.extend([pcedim]*uncdims)
    ypcedims = tuple(ypcedims)
    ytens = np.zeros(ypcedims)

    abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)
    # abscarray, weightsarray = np.array(abscissae), np.array(weights)

    nulist = [basenu]*5
    nuarray = np.array(nulist)
    for idxtuple in itertools.product(np.arange(pcedim), repeat=uncdims):
        idxarray = np.array(idxtuple)
        nuarray[:uncdims] = abscissae[idxarray]
        pceylist.append(get_output(nuarray.tolist(), plotfignum=None))

    yrslttns = np.array(pceylist).reshape(ypcedims)
    yrslttns = tnsrtrnsps(yrslttns)
    exypce = 0
    # for kk, cw in enumerate(weights):
    for idxtuple in itertools.product(np.arange(pcedim), repeat=uncdims):
        idxarray = np.array(idxtuple)
        cw = (weights[idxarray]).prod()
        exypce += cw*yrslttns[idxtuple]

    print('pcedim={0:2.0f}, exypce={1}'.
          format(pcedim, 1./((varib-varia)**uncdims)*exypce))
