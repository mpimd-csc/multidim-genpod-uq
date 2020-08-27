import numpy as np
from mc_pce_gp import simit

problem = 'cylinder'
plotplease = False
meshlevellist = np.arange(5, 12)

dofslist, ylist = [], []
for meshlevel in meshlevellist:
    dofs, outpt = simit(problem=problem, meshlevel=meshlevel,
                        nulb=4e-4, nuub=4e-4,
                        plotplease=plotplease, onlymeshtest=True)
    dofslist.append(dofs)
    ylist.append(outpt)

for k, ml in enumerate(meshlevellist):
    print('Mesh:{0} | dofs:{1} | y:{2}'.format(ml, dofslist[k], ylist[k]))
