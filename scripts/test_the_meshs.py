import numpy as np
from mc_pce_gp import simit
import logging
from rich.logging import RichHandler

logging.basicConfig(level=logging.INFO, handlers=[RichHandler()],
                    format='%(message)s', datefmt="[%X]",)

problem = 'cylinder'
plotplease = False
meshlevellist = np.arange(5, 12)
meshlevellist = np.arange(12, 13)
meshlevellist = np.arange(4, 8)
distribution = 'beta-2-5'
multiproc = 4
fullsweep = True

simpars = dict(problem=problem, multiproc=multiproc,
               distribution=distribution,
               omtp_dict=dict(fullsweep=fullsweep, pcedim=2),
               nulb=5e-4, nuub=1e-3, plotplease=plotplease, onlymeshtest=True)

dofslist, ylist = [], []
for meshlevel in meshlevellist:
    simpars.update(dict(meshlevel=meshlevel))
    if fullsweep:
        dofs, outpt, absc = simit(**simpars)
        ylist.append(outpt.reshape(-1))
    else:
        dofs, outpt = simit(**simpars)
    dofslist.append(dofs)

if fullsweep:
    np.set_printoptions(precision=4)
    for kkk, ml in enumerate(meshlevellist[1:]):
        valdiff = ylist[kkk+1] - ylist[kkk]
        print(f'Mesh:{ml} | {valdiff}')
else:
    for k, ml in enumerate(meshlevellist):
        b = [f'{x:.5f}' for x in ylist[k]]
        print('Mesh:{0} | dofs:{1} | y:{2}'.format(ml, dofslist[k], b))
