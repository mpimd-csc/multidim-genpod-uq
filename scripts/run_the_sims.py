import getopt
import sys
import numpy as np
import re

from mc_pce_gp import simit
from dolfin_navier_scipy.data_output_utils import Timer

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
pceplease = True
# plotplease = True
pcepod = True
# mcpod = True

basisfrom = 'mc'
basisfrom = 'pce'
problem = 'cylinder'
meshlevel = 10
nulb = 4e-4
nuub = 6e-4
nprocs = 2
timings = 1

plotpcepoddiff = False
pcepoddiffdim = 9

mcxpy, pcexpy = None, None


options, rest = getopt.getopt(sys.argv[1:], '',
                              ['mesh=',
                               'mc=',
                               'mcruns=',
                               'pce=',
                               'pcedims=',
                               'poddims=',
                               'podbase=',
                               'pcepod=',
                               'mcpod=',
                               'pcesnapdim=',
                               'mcsnap=',
                               'varinu=',
                               'nprocs=',
                               'mcxpy=',
                               'pcexpy=',
                               'timings=',
                               'plotpcepoddifddim='
                               ])

for opt, arg in options:
    if opt == '--mesh':
        meshlevel = np.int(arg)
    elif opt == '--pcepod':
        pcepod = np.bool(np.int(arg))
    elif opt == '--mcpod':
        mcpod = np.bool(np.int(arg))
    elif opt == '--mc':
        mcplease = np.bool(np.int(arg))
    elif opt == '--mcruns':
        mcruns = np.int(arg)
        if mcruns < 10:
            raise UserWarning('minimal number for mcruns is 10')
    elif opt == '--podbase':
        basisfrom = np.str(arg)
    elif opt == '--pce':
        pceplease = np.bool(np.int(arg))
    elif opt == '--pcedims':
        nmstrl = re.findall('\\d+', arg)
        pcedimlist = [np.int(xstr) for xstr in nmstrl]
    elif opt == '--poddims':
        poddstrl = re.findall('\\d+', arg)
        poddimlist = [np.int(xstr) for xstr in poddstrl]
    elif opt == '--pcesnapdim':
        pcesnapdim = np.int(arg)
    elif opt == '--mcsnap':
        mcsnap = np.int(arg)
    elif opt == '--varinu':
        nuabstrl = re.findall('\\d+', arg)
        nuabl = [np.int(xstr) for xstr in nuabstrl]
        nulb = nuabl[0]*10**(-nuabl[2])
        nuub = nuabl[1]*10**(-nuabl[2])
    elif opt == '--nprocs':
        nprocs = np.int(arg)
    elif opt == '--mcxpy':
        mcxpy = np.float(arg)
    elif opt == '--pcexpy':
        pcexpy = np.float(arg)
    elif opt == '--timings':
        timings = np.int(arg)
    elif opt == '--plotpcepoddifddim':
        pcepoddiffdim = np.int(arg)
        if pcepoddiffdim > 0:
            plotpcepoddiff = True

infostring = ('meshlevel      = {0}'.format(meshlevel) +
              '\npce            = {0}'.format(pceplease) +
              '\nmc             = {0}'.format(mcplease) +
              '\nmcpod          = {0}'.format(mcpod) +
              '\npcepod         = {0}'.format(pcepod)
              )

if mcplease:
    infostring = (infostring +
                  '\nmcruns         = {0}'.format(mcruns))

if pceplease or pcepod:
    infostring = (infostring +
                  '\npcedims        = {0}'.format(pcedimlist))

if mcpod or pcepod:
    infostring = (infostring +
                  '\npoddimlist     = {0}'.format(poddimlist) +
                  '\nbasisfrom      = {0}'.format(basisfrom))
    if basisfrom == 'mc':
        infostring = (infostring +
                      '\nmc snapshots   = {0}'.format(mcsnap) +
                      '\nred mc runs    = {0}'.format(mcruns))
    if basisfrom == 'pce':
        infostring = (infostring +
                      '\ntrain pce dim  = {0}'.format(pcesnapdim))
if nprocs > 1:
    infostring = (infostring +
                  '\nnprocs         = {0}'.format(nprocs))


print('******************')
print(infostring)
print('******************')

with Timer():
    simit(mcruns=mcruns, pcedimlist=pcedimlist,
          problem=problem, meshlevel=meshlevel,
          plotplease=plotplease, basisfrom=basisfrom,
          mcxpy=mcxpy, pcexpy=pcexpy, redmcruns=15000,
          mcsnap=mcsnap, pcesnapdim=pcesnapdim, poddimlist=poddimlist,
          multiproc=nprocs, timings=timings,
          # basenu=basenu, varia=varia, varib=varib,
          plotpcepoddiff=plotpcepoddiff, pcepoddiffdim=pcepoddiffdim,
          nulb=nulb, nuub=nuub,
          mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
