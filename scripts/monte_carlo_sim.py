import numpy as np
import matplotlib.pyplot as plt

from circle_subsec import get_problem

get_linop, get_sol, get_output, problemfems = get_problem()

basenu = 1e-3

Nruns = 100


def realizenulist(poslist):
    nulist = [basenu]*5
    for k in poslist:
        nulist[k] = basenu + ((1e-4*np.random.rand(1)).flatten())[0]
    return nulist

estxnu, estxy = 0, 0
ylist, xnulist = [], []
for mck in range(Nruns):
    nulist = realizenulist([4])
    estxnu += nulist[4]
    cury = get_output(nulist, plotfignum=None)
    ylist.append(cury)
    xnulist.append(nulist[4])
    estxy += cury

estxnu = 1./Nruns*estxnu
estxy = 1./Nruns*estxy

print('estxnu={0}, estxy={1}'.format(estxnu, estxy))

nulist = [basenu]*5
nulist[4] = estxnu
cury = get_output(nulist, plotfignum=None)
print('y(estxnu)={0}'.format(cury))

plt.figure(89)
plt.plot(ylist, '.')
plt.figure(98)
plt.plot(xnulist, '.')
plt.show()
