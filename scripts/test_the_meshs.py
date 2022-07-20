import numpy as np
import re
from mc_pce_gp import simit
import logging
import argparse
from rich.logging import RichHandler

logging.basicConfig(level=logging.INFO, handlers=[RichHandler()],
                    format='%(message)s', datefmt="[%X]",)

problem = 'cylinder'
plotplease = False
fullsweep = False
smlmesh = 7
lrgmesh = 8
if fullsweep:
    distribution = 'beta-2-5'
    distribution = 'uniform'
else:
    distribution = None  # not needed
multiproc = 4

nulb = 5e-4
nuub = 1e-3

prsr = argparse.ArgumentParser()
prsr.add_argument("--distribution", type=str, help="type of distribution",
                  default=distribution)
prsr.add_argument("--pcetestdim", type=int,
                  help="dimensions PCE for the mesh test", default=2)
prsr.add_argument("--varinu", type=str, help="range of the nu's", default=None)
prsr.add_argument("--nprocs", type=int,
                  help="number of parallel threads", default=multiproc)
prsr.add_argument("--smlmesh", type=int,
                  help="lowest dimension of mesh to test", default=smlmesh)
prsr.add_argument("--lrgmesh", type=int,
                  help="lowest dimension of mesh to test", default=lrgmesh)
args = prsr.parse_args()
logging.info(args)
meshlevellist = np.arange(args.smlmesh, args.lrgmesh+1)

if args.varinu is not None:
    nuabstrl = re.findall('\\d+', args.varinu)
    nuabl = [np.int(xstr) for xstr in nuabstrl]
    nulb = nuabl[0]*10**(-nuabl[2])
    nuub = nuabl[1]*10**(-nuabl[2])
else:
    pass

simpars = dict(problem=problem, multiproc=args.nprocs,
               distribution=args.distribution,
               omtp_dict=dict(fullsweep=fullsweep, pcedim=args.pcetestdim),
               nulb=nulb, nuub=nuub, plotplease=plotplease, onlymeshtest=True)

dofslist, ylist = [], []
for meshlevel in meshlevellist:
    simpars.update(dict(meshlevel=meshlevel))
    if fullsweep:
        dofs, outpt, absc = simit(**simpars)
        ylist.append(outpt.reshape(-1))
    else:
        dofs, outpt = simit(**simpars)
        ylist.append(outpt)
    dofslist.append(dofs)

if fullsweep:
    np.set_printoptions(precision=4)
    for kkk, ml in enumerate(meshlevellist[1:]):
        valdiff = ylist[kkk+1] - ylist[kkk]
        print(f'mshd:{ml}-{meshlevellist[kkk]} | {valdiff}')
else:
    for k, ml in enumerate(meshlevellist):
        b = [f'{x:.5f}' for x in ylist[k]]
        print('Mesh:{0} | dofs:{1} | y:{2}'.format(ml, dofslist[k], b))
