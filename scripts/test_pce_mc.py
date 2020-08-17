import numpy as np

import gen_pod_uq.mc_pce_utils as mpu

nulb = 4e-4
nuub = 6e-4

nprocs = 1
uncdims = 1

pcedim = 6
mcruns = 1000000

varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
varinulist = varinu.tolist()


def get_output(alpha):
    alpharr = np.array(alpha)
    return 1./alpharr


mcout, mcxpy, expvnu = mpu.run_mc_sim(varinulist, get_output,
                                      verbose=True,
                                      multiproc=nprocs)

abscissae, weights, compexpv = mpu.\
    setup_pce(distribution='uniform',
              distrpars=dict(a=nulb, b=nuub),
              pcedim=pcedim, uncdims=uncdims)
ysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
                                     uncdims=uncdims,
                                     multiproc=nprocs,
                                     abscissae=abscissae)
pcexpy = compexpv(ysoltens)
trvl = (np.log(nuub) - np.log(nulb))/(nuub - nulb)
print('PCE({0}): E(y): {1}'.format(pcedim, pcexpy))
print('MC({0}): E(y): {1}'.format(mcruns, mcxpy))
print('true value: {0}'.format(trvl))
