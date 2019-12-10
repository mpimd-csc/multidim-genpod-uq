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
meshlevel = 1

options, rest = getopt.getopt(sys.argv[1:], '',
                              ['mesh='
                               ])
for opt, arg in options:
    if opt == '--mesh':
        meshlevel = np.int(arg)


simit(mcruns=mcruns, pcedimlist=pcedimlist,
      problem=problem, meshlevel=meshlevel,
      plotplease=plotplease, basisfrom=basisfrom,
      mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
