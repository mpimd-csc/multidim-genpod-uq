import getopt
import sys
import numpy as np

from mc_pce_gp import simit

mcruns = 10000  # 200
pcedimlist = [2]  # , 5]  # , 3, 4, 5]  # , 7]
mcplease = False
pceplease = False
plotplease = False
mcpod = False
pcepod = False
# ## make it come true
# mcplease = True
# pceplease = True
plotplease = True
# pcepod = True
# mcpod = True
basisfrom = 'mc'
basisfrom = 'pce'
problem = 'cylinder'
meshlevel = 6

options, rest = getopt.getopt(sys.argv[1:], '',
                              ['mesh=',
                               'mc=',
                               'pce=',
                               'podbas=',
                               'pcepod=',
                               'mcpod='
                               ])

for opt, arg in options:
    if opt == '--mesh':
        meshlevel = np.int(arg)
    elif opt == '--pcepod':
        pcepod = np.bool(np.int(arg))
    elif opt == '--mcpod':
        mcpod = np.bool(np.int(arg))
    elif opt == '--mc':
        mcruns = np.int(arg)
        if mcruns >= 10:
            mcplease = True
        else:
            print('minimal number for mcruns is 10')
            mcplease = False
    elif opt == '--podbas':
        basisfrom = np.str(arg)

infostring = ('meshlevel      = {0}'.format(meshlevel) +
              '\nbasisfrom      = {0}'.format(basisfrom) +
              '\npce            = {0}'.format(pceplease) +
              '\nmc             = {0}'.format(mcplease)
              )

print(infostring)


simit(mcruns=mcruns, pcedimlist=pcedimlist,
      problem=problem, meshlevel=meshlevel,
      plotplease=plotplease, basisfrom=basisfrom,
      mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
