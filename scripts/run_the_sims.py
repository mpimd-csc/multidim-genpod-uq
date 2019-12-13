import getopt
import sys
import numpy as np

from mc_pce_gp import simit

mcruns = 10000  # 200
pcedimlist = [3, 4]  # , 5]  # , 3, 4, 5]  # , 7]
mcplease = False
pceplease = False
plotplease = False
mcpod = False
pcepod = False
# ## make it come true
# mcplease = True
# pceplease = True
# plotplease = True
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
                               'podbase=',
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
    elif opt == '--podbase':
        basisfrom = np.str(arg)
    elif opt == '--pce':
        pceplease = np.bool(np.int(arg))

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

print('******************')
print(infostring)
print('******************')


simit(mcruns=mcruns, pcedimlist=pcedimlist,
      problem=problem, meshlevel=meshlevel,
      plotplease=plotplease, basisfrom=basisfrom,
      mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
