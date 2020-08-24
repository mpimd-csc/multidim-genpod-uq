import numpy as np

import gen_pod_uq.mc_pce_utils as mpu

nulb = 3e-4
nuub = 7e-4
uncdims = 1
uncdims = 2
cplpar = 1e-4

nprocs = 1

pcedims = [3, 4, 5, 6]
mcrunslist = [10**x for x in range(4, 7)]

mcsims = 5

if uncdims == 1:
    def get_output(alpha):
        alpharr = np.array(alpha)
        return 1./alpharr

elif uncdims == 2:
    def get_output(alpha):
        alpharr = np.array(alpha)
        aone, atwo = alpharr[0], alpharr[1]
        return 1./(aone*atwo - cplpar**2) * (aone + atwo - 2*cplpar)

# Used Mathematica for these values (see the code below)
if uncdims == 1:
    trvl = 2118.24465097
    vtrvl = 274944.360550
elif uncdims == 2:
    trvl = 3504.22709343
    vtrvl = 261037.034256

# ## PCE approximation
pcexps, pcevrs = [], []
for pcedim in pcedims:
    abscissae, weights, compexpv, compvrnc = mpu.\
        setup_pce(distribution='uniform',
                  distrpars=dict(a=nulb, b=nuub),
                  pcedim=pcedim, uncdims=uncdims)
    ysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
                                         uncdims=uncdims,
                                         multiproc=nprocs,
                                         abscissae=abscissae)
    pcexpy = compexpv(ysoltens)
    pcvrnc = compvrnc(ysoltens, pcexpy)
    pcexps.append(pcexpy)
    pcevrs.append(pcvrnc)

mcexps, mcevrs = [], []
fmnusqrd = mpu.solfunc_to_variance(get_output, trvl)
for mcruns in mcrunslist:
    mcxpvs = []
    mcvrcs = []
    for nmcsim in range(mcsims):
        varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
        varinulist = varinu.tolist()
        mcvarout, mcvar, expvnu = mpu.run_mc_sim(varinulist, fmnusqrd,
                                                 multiproc=nprocs)
        mcout, mcxpy, expvnu = mpu.run_mc_sim(varinulist, get_output,
                                              multiproc=nprocs)
        mcxpvs.append(mcxpy)
        mcvrcs.append(mcvar)

    mcexps.append(np.median(np.array(mcxpvs)))
    mcevrs.append(np.median(np.array(mcvrcs)))

pcexps = np.array(pcexps)
pcexpserrs = (pcexps - trvl)/trvl
pcevrs = np.array(pcevrs)
pcvrserrs = (pcevrs - vtrvl)/vtrvl

mcexps = np.array(mcexps)
mcexpserrs = (mcexps - trvl)/trvl
mcevrs = np.array(mcevrs)
mcvrserrs = (mcevrs - vtrvl)/vtrvl

print('pcedims:', pcedims)
print('rel err expv:', pcexpserrs)
print('rel err vrnc:', pcvrserrs)

print('mcruns:', mcrunslist)
print('rel err expv:', mcexpserrs)
print('rel err vrnc:', mcvrserrs)

# print('true value: E(y): {0}'.format(trvl))
# print('true value: V(y): {0}'.format(vtrvl))
# print('PCE({0}): E(y): {1}'.format(pcedim, pcexpy))
# print('PCE({0}): V(y): {1}'.format(pcedim, pcvrnc))
# print('MC({0}): E(y): {1}'.format(mcruns, mcxpy))
# print('MC({0}): V(y): {1}'.format(mcruns, mcvar))
# print('MC E(Y): {0} -- median out of {1}'.format(medexpv, mcsims))
# print('MC V(Y): {0} -- median out of {1}'.format(medvrnc, mcsims))

'''
Mathematica code for the true values

eps = 10^-4
ao = 3*10^-4
bo = 7*10^-4
at = ao
bt = bo
sfo = 1/(bo - ao)
fo = 1/x
sft = 1/((bo - ao)*(bt - at))
ft = 1/(x*y - eps^2)*(x + y - 2*eps)
exo = Integrate[sfo*fo, {x, ao, bo}]
vxo = Integrate[sfo*fo^2, {x, ao, bo}]
ext = Integrate[sft*ft, {x, ao, bo}, {y, at, bt}]
vxt = Integrate[sft*ft^2, {x, ao, bo}, {y, at, bt}]
N[exo, 12]
N[vxo - exo^2, 12]
N[ext, 12]
N[vxt - ext^2, 12]
'''
