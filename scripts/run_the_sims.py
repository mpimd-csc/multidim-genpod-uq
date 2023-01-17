# import getopt
import argparse
# import numpy as np
import re

from mc_pce_gp import simit
from dolfin_navier_scipy.data_output_utils import Timer

import logging
from rich.logging import RichHandler

logging.basicConfig(level=logging.INFO, handlers=[RichHandler()],
                    format='%(message)s',
                    datefmt="[%X]",
                    )

mcruns = 1000  # 200
pcedimlist = [3, 4, 5]  # , 3, 4, 5]  # , 7]
pcesnapdim = 3
mcsnap = 3**5*2
poddimlist = [5, 10, 20]

mcplease = False
pceplease = False
plotplease = False
mcpod = False
pcepod = False
# ## make it come true
# mcplease = True
# pceplease = True
# plotplease = True
pcepod = True
# mcpod = True

basisfrom = 'mc'
basisfrom = 'pce'
problem = 'cylinder'
meshlevel = 6
nulb = 3e-4
nuub = 7e-4
distribution = 'beta-2-5'
distribution = 'uniform'
nprocs = 4
timings = 1

plotpcepoddiff = False
pcepoddiffdim = 9

mcxpy, pcexpy = None, None

prsr = argparse.ArgumentParser()
prsr.add_argument("--mesh", type=int, help="mesh level", default=meshlevel)
prsr.add_argument("--pcepod", type=lambda x: bool(int(x)),
                  help="do POD PCE", default=pcepod)
prsr.add_argument("--distribution", type=str, help="type of distribution",
                  default=distribution)
prsr.add_argument("--mcpod", type=lambda x: bool(int(x)),
                  help="do POD MC", default=mcpod)
prsr.add_argument("--mc", type=lambda x: bool(int(x)),
                  help="do FOM MC", default=mcplease)
prsr.add_argument("--mcruns", type=int, help="#MC samples", default=mcruns)
prsr.add_argument("--rombase", type=str,
                  help="what basis for the ROM", default=basisfrom)
prsr.add_argument("--pce", type=lambda x: bool(int(x)),
                  help="do a FOM PCE", default=pceplease)
prsr.add_argument("--pcedims", type=str,
                  help="dimensions of the FOM/ROM PCE", default=None)
prsr.add_argument("--poddims", type=str, help="dimensions ROM", default=None)
prsr.add_argument("--pcesnapdim", type=int,
                  help="dimensions of training PCE", default=pcesnapdim)
prsr.add_argument("--mcsnap", type=int,
                  help="dimensions of training MC", default=mcsnap)
prsr.add_argument("--rbsnap", type=int,
                  help="dimensions of training set for RB", default=None)
prsr.add_argument("--varinu", type=str, help="range of the nu's", default=None)
prsr.add_argument("--nprocs", type=int,
                  help="number of parallel threads", default=nprocs)
prsr.add_argument("--timings", type=int,
                  help="number of runs for timing", default=timings)

args = prsr.parse_args()
logging.info(args)

if args.pcedims is not None:
    nmstrl = re.findall('\\d+', args.pcedims)
    pcedimlist = [int(xstr) for xstr in nmstrl]
else:
    pass

if args.poddims is not None:
    poddstrl = re.findall('\\d+', args.poddims)
    poddimlist = [int(xstr) for xstr in poddstrl]
else:
    pass

if args.varinu is not None:
    nuabstrl = re.findall('\\d+', args.varinu)
    nuabl = [int(xstr) for xstr in nuabstrl]
    nulb = nuabl[0]*10**(-nuabl[2])
    nuub = nuabl[1]*10**(-nuabl[2])
else:
    pass

if args.mcruns < 10:
    raise UserWarning('minimal number for mcruns is 10')
else:
    pass

infostring = ('meshlevel      = {0}'.format(args.mesh) +
              '\ndistribution   = ' + args.distribution +
              '\npce            = {0}'.format(args.pce) +
              '\nmc             = {0}'.format(args.mc) +
              '\nmcpod          = {0}'.format(args.mcpod) +
              '\npcepod         = {0}'.format(args.pcepod) +
              f'\nbasisfrom      = {args.rombase}'
              )

if args.mc:
    infostring = (infostring +
                  '\nmcruns         = {0}'.format(args.mcruns))

if args.pce or args.pcepod:
    infostring = (infostring +
                  '\npcedims        = {0}'.format(pcedimlist))

if args.mcpod or args.pcepod:
    infostring = (infostring +
                  '\npoddimlist     = {0}'.format(poddimlist) +
                  '\nbasisfrom      = {0}'.format(args.rombase))
    if args.rombase == 'mc':
        infostring = (infostring +
                      '\nmc snapshots   = {0}'.format(args.mcsnap) +
                      '\nred mc runs    = {0}'.format(args.mcruns))
    elif args.rombase == 'pce':
        infostring = (infostring +
                      '\ntrain pce dim  = {0}'.format(args.pcesnapdim))
    elif args.rombase == 'rb':
        infostring = (infostring +
                      '\ntrain rb dim   = {0}'.format(args.rbsnap))
    else:
        pass
else:
    pass

if nprocs > 1:
    infostring = (infostring +
                  '\nnprocs         = {0}'.format(nprocs))
else:
    pass

logging.info('******************')
logging.info(infostring)
logging.info('******************')

with Timer():
    simit(mcruns=args.mcruns, pcedimlist=pcedimlist,
          problem=problem, meshlevel=args.mesh,
          plotplease=plotplease, basisfrom=args.rombase,
          # mcxpy=args.mcxpy, pcexpy=args.pcexpy,
          # redmcruns=15000,
          mcsnap=args.mcsnap, trainpcedim=args.pcesnapdim,
          poddimlist=poddimlist,
          rbparams=dict(samplemethod='random', nsample=args.rbsnap),
          multiproc=args.nprocs, timings=args.timings,
          # plotpcepoddiff=args.plotpcepoddiff, pcepoddiffdim=pcepoddiffdim,
          nulb=nulb, nuub=nuub,
          distribution=args.distribution,
          mcplease=args.mc, pceplease=args.pce,
          mcpod=args.mcpod, pcepod=args.pcepod)
