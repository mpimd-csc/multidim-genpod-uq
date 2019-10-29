import numpy as np
import matplotlib.pyplot as plt

import spacetime_galerkin_pod.chaos_expansion_utils as ceu

from circle_subsec import get_problem

get_linop, get_sol, get_output, problemfems = get_problem()

uncdims = 2

basenu = 1e-3
varia = 0.
varib = 1e-4

mcits, mcruns = 5, 5  # 200
pcedimlist = [3, 5, 8, 12, 17]

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

cury = get_output(nulist, plotfignum=None)
print('y(estxnu)={0}'.format(cury))

if plotplease:
    plt.figure(89)
    plt.plot(ylist, '.')
    # plt.figure(98)
    # plt.plot(xnulist, '.')
    plt.show()


# ## CHAP Polynomial Chaos Expansion

nua, nub = basenu+varia, basenu+varib

for pcedim in pcedimlist:
    abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)

    nulist = [basenu]*5
    ylist = []
    for cnu in abscissae:
        nulist[0] = cnu
        ylist.append(get_output(nulist, plotfignum=None))

    exypce = 0
    for kk, cw in enumerate(weights):
        exypce += cw*ylist[kk]

    print('pcedim={0:2.0f}, exypce={1}'.format(pcedim,
                                               1./(varib-varia)*exypce))
