import itertools

import numpy as np

import spacetime_galerkin_pod.chaos_expansion_utils as ceu


def run_mc_sim(parlist, solfunc, chunks=10,
               comp_para_ev=True, verbose=False, ret_ev=True):
    expvpara = None
    if comp_para_ev:
        expvpara = np.average(np.array(parlist), axis=0)

    ylist = []
    if verbose:
        nmc = len(parlist)
        mcx = np.int(nmc/chunks)
        for mcit in range(chunks-1):
            for cpar in parlist[mcit*mcx:(mcit+1)*mcx]:
                cy = (solfunc(cpar)).flatten()
                ylist.append(cy)
            estxy = np.average(np.array(ylist), axis=0)[0]
            print('mc:{0}/{1}: estxy[0]={2}'.
                  format((mcit+1)*mcx, nmc, estxy))
        for cpar in parlist[(mcit+1)*mcx:]:
            cy = (solfunc(cpar)).flatten()
            ylist.append(cy)
        estxy = np.average(np.array(ylist), axis=0)[0]
        print('mc:{0}/{1}: estxy[0]={2}'.format(nmc, nmc, estxy))

        return ylist, np.average(np.array(ylist), axis=0), expvpara

    else:
        for cpar in parlist:
            ylist.append(solfunc(cpar))
        return ylist, np.average(np.array(ylist), axis=0), expvpara


def run_pce_sim_separable(solfunc=None, uncdims=None, abscissae=None):
    """ pce simulation for all PCE dimensions being the same
    """
    # compute the sols
    pceylist = []
    for absctpl in itertools.product(abscissae, repeat=uncdims):
        pceylist.append((solfunc(absctpl)).flatten())

    ypcedims = [pceylist[0].size]
    ypcedims.extend([len(abscissae)]*uncdims)
    ypcedims = tuple(ypcedims)
    # arrange it in the tensor
    # first dimension is the state, then comes the uncertainty
    yrslttns = (np.array(pceylist).T).reshape(ypcedims)
    return yrslttns


def setup_pce(distribution='uniform', distrpars={}, pcedim=None, uncdims=None):
    # compute the expected value
    if distribution == 'uniform':
        abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, **distrpars)

    scalefac = (1./weights.sum())**uncdims

    def comp_expv(ytens):
        ydim = ytens.shape[0]
        expv = 0
        maty = ytens.reshape((ydim, -1))
    # for kk, cw in enumerate(weights):
        for idx, wtpl in enumerate(itertools.product(weights, repeat=uncdims)):
            cw = (np.array(wtpl)).prod()
            expv += cw*maty[:, idx]
        return scalefac*expv

    return abscissae, weights, comp_expv
