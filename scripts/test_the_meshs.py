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
distribution = 'beta-2-5'

dofslist, ylist = [], []
for meshlevel in meshlevellist:
    dofs, outpt = simit(problem=problem, meshlevel=meshlevel,
                        distribution=distribution,
                        nulb=3e-4, nuub=7e-4,
                        plotplease=plotplease, onlymeshtest=True)
    dofslist.append(dofs)
    ylist.append(outpt)

for k, ml in enumerate(meshlevellist):
    b = [f'{x:.5f}' for x in ylist[k]]
    print('Mesh:{0} | dofs:{1} | y:{2}'.format(ml, dofslist[k], b))
