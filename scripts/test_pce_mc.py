import numpy as np

import gen_pod_uq.mc_pce_utils as mpu

from mat_lib_plots import conv_plot_utils as cpu

nulb = 3e-4
nuub = 7e-4
cplpar = 1e-4

nprocs = 1


def simit(uncdims=1, mcsims=5, pcedims=[3],
          mcrunslist=[100000]):

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
        pcexpysqrd = compexpv(np.square(ysoltens))
        pcexps.append(pcexpy)
        pcevrs.append(pcexpysqrd-pcexpy**2)

    mcexps, mcevrs = [], []
    for mcruns in mcrunslist:
        mcxpvs = []
        mcvrcs = []

        for nmcsim in range(mcsims):
            varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
            varinulist = varinu.tolist()
            mcout, mcxpy, expvnu = mpu.run_mc_sim(varinulist, get_output,
                                                  multiproc=nprocs)
            mcxpvs.append(mcxpy)
        mcexps.append(np.median(np.array(mcxpvs)))

        for nmcsim in range(mcsims):
            varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
            varinulist = varinu.tolist()
            fmnusqrd = mpu.solfunc_to_variance(get_output, mcexps[-1])
            mcvarout, mcvar, expvnu = mpu.run_mc_sim(varinulist, fmnusqrd,
                                                     multiproc=nprocs)
            mcvrcs.append(mcvar)

        mcevrs.append(np.median(np.array(mcvrcs)))

    pcexps = np.array(pcexps)
    pcexpserrs = (pcexps - trvl)/trvl
    pcevrs = np.array(pcevrs)
    pcvrserrs = (pcevrs - vtrvl)/vtrvl

    mcexps = np.array(mcexps)
    mcexpserrs = (mcexps - trvl)/trvl
    mcevrs = np.array(mcevrs)
    mcvrserrs = (mcevrs - vtrvl)/vtrvl

    pcenms = ['\\pcex {0}'.format(x) for x in pcedims]
    mcnms = ['\\mcx {0}'.format(x) for x in mcrunslist]

    print('\n ****************')
    print('### {0}D Problem ###'.format(uncdims))
    print(' ****************')
    # cpu.print_nparray_tex(np.array(pcexpserrs), name=pcenms, fstr='.2e')
    print('\n PCE - expv - rel err:')
    cpu.print_nparray_tex(np.array(pcexpserrs).reshape((1, len(pcedims))),
                          name=pcenms, fstr='.2e')
    print('\n PCE - vrnc - rel err:')
    cpu.print_nparray_tex(np.array(pcvrserrs).reshape((1, len(pcedims))),
                          name=pcenms, fstr='.2e')
    print('\n MC - expv - rel err:')
    cpu.print_nparray_tex(np.array(mcexpserrs).reshape((1, len(mcrunslist))),
                          name=mcnms, fstr='.2e')
    print('median out of {0} runs'.format(mcsims))

    print('\n MC - vrnc - rel err:')
    cpu.print_nparray_tex(np.array(mcvrserrs).reshape((1, len(mcrunslist))),
                          name=mcnms, fstr='.2e')
    print('median out of {0} runs'.format(mcsims))


if __name__ == '__main__':
    pcedims = [3, 4, 5, 6]
    mcrunslist = [10**x for x in range(4, 7)]
    mcsims = 15

    np.random.seed(1)
    simit(uncdims=1, pcedims=pcedims, mcrunslist=mcrunslist, mcsims=mcsims)

    np.random.seed(1)
    simit(uncdims=2, pcedims=pcedims, mcrunslist=mcrunslist, mcsims=mcsims)


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
