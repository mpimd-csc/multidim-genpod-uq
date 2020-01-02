import getopt
import sys
import numpy as np

from mc_pce_gp import simit

mcruns = 1000  # 200
pcedimlist = [3, 4]  # , 5]  # , 3, 4, 5]  # , 7]
pcesnapdim = 3
mcsnap = 3**5*2
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
mcpod = True
basisfrom = 'mc'
basisfrom = 'pce'
problem = 'cylinder'
meshlevel = 1

if meshlevel == 6:
    mcxpy = 0.7234999474635652  # from 15000 mc runs, see editha-logs
    pcexpy = 0.72347945  # PCE(5) see editha-logs
elif meshlevel == 1:
    mcxpy = 0.660508808729938  # from 1000 mc runs
    pcexpy = 0.66298159  # PCE(4)
else:
    mcxpy, pcexpy = None, None

options, rest = getopt.getopt(sys.argv[1:], '',
                              ['mesh=',
                               'mc=',
                               'mcruns=',
                               'pce=',
                               'podbase=',
                               'pcepod=',
                               'mcpod=',
                               'pcesnapdim=',
                               'mcsnap='
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
    elif opt == '--pcesnapdim':
        pcesnapdim = np.int(arg)
    elif opt == '--mcsnap':
        mcsnap = np.int(arg)

infostring = ('meshlevel      = {0}'.format(meshlevel) +
              '\nbasisfrom      = {0}'.format(basisfrom) +
              '\npce            = {0}'.format(pceplease) +
              '\nmc             = {0}'.format(mcplease) +
              '\nmcpod          = {0}'.format(mcpod) +
              '\npcepod         = {0}'.format(pcepod)
              )

if mcplease:
    infostring = (infostring +
                  '\nmcruns         = {0}'.format(mcruns))

if basisfrom == 'mc':
    infostring = (infostring +
                  '\nmc snapshots   = {0}'.format(mcsnap) +
                  '\nred mc runs    = {0}'.format(mcruns))

if basisfrom == 'pce':
    infostring = (infostring +
                  '\ntrain pce dim  = {0}'.format(pcesnapdim))

print('******************')
print(infostring)
print('******************')


simit(mcruns=mcruns, pcedimlist=pcedimlist,
      problem=problem, meshlevel=meshlevel,
      plotplease=plotplease, basisfrom=basisfrom,
      mcxpy=mcxpy, pcexpy=pcexpy, redmcruns=15000,
      mcsnap=mcsnap, pcesnapdim=pcesnapdim,
      mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
