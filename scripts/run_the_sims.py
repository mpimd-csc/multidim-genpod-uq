from mc_pce_gp import simit

mcruns = 10  # 200
pcedimlist = [3]  # , 5]  # , 3, 4, 5]  # , 7]
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
meshlevel = 1

simit(mcruns=mcruns, pcedimlist=pcedimlist,
      problem=problem, meshlevel=meshlevel,
      plotplease=plotplease, basisfrom=basisfrom,
      mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
