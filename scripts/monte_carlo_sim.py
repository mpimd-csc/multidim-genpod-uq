import numpy as np
import matplotlib.pyplot as plt

import spacetime_galerkin_pod.chaos_expansion_utils as ceu

from circle_subsec import get_problem

get_linop, get_sol, get_output, problemfems = get_problem()

basenu = 1e-3
varia = 0.
varib = 1e-4

mcits, mcruns = 25, 1000  # 200
pcedimlist = [3, 5, 8, 12, 17]

plotplease = False

# ## CHAP Monte Carlo


def realizenulist(poslist):
    nulist = [basenu]*5
    for k in poslist:
        nulist[k] = basenu + (((varib-varia)*np.random.rand(1)).flatten())[0]
    return nulist

# estxnu, estxy = 0, 0
ylist, xnulist = [], []
for mitk in range(mcits):
    for mck in range(mcruns):
        nulist = realizenulist([4])
        # estxnu += nulist[4]
        cury = get_output(nulist, plotfignum=None)
        ylist.append(cury)
        xnulist.append(nulist[4])
        # estxy += cury
    estxnu = np.average(np.array(nulist))
    estxy = np.average(np.array(ylist))

    print('mc:{0}/{1}: estxy={2}'.format((mitk+1)*mcruns, mcits*mcruns, estxy))

nulist = [basenu]*5
nulist[4] = basenu+.5*(varib-varia)  # estxnu
cury = get_output(nulist, plotfignum=None)
print('y(estxnu)={0}'.format(cury))

if plotplease:
    plt.figure(89)
    plt.plot(ylist, '.')
    plt.figure(98)
    plt.plot(xnulist, '.')
    plt.show()


# ## CHAP Polynomial Chaos Expansion

nua, nub = basenu+varia, basenu+varib

for pcedim in pcedimlist:
    abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)

    nulist = [basenu]*5
    ylist = []
    for cnu in abscissae:
        nulist[4] = cnu
        ylist.append(get_output(nulist, plotfignum=None))

    exypce = 0
    for kk, cw in enumerate(weights):
        exypce += cw*ylist[kk]

    print('pcedim={0:2.0f}, exypce={1}'.format(pcedim,
                                               1./(varib-varia)*exypce))
