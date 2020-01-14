from itertools import product, islice
from multiprocessing import Queue, Process

import numpy as np

import spacetime_galerkin_pod.chaos_expansion_utils as ceu


def run_mc_sim(parlist, solfunc, chunks=10, multiproc=0,
               comp_para_ev=True, verbose=False, ret_ev=True):
    expvpara = None
    if comp_para_ev:
        expvpara = np.average(np.array(parlist), axis=0)

    def _compallsols(paralist):
        locylist = []
        for cpara in paralist:
            locylist.append(solfunc(cpara))
        return locylist

    def _formproc(paralist, pqueue):
        ylist = _compallsols(paralist)
        pqueue.put(ylist)

    nmc = len(parlist)
    if multiproc > 1:
        pqueue = Queue()
        mcx = nmc/multiproc

        itschunks = []
        for k in range(multiproc-1):
            itschunks.append(parlist[np.int(np.floor(k*mcx)):
                                     np.int(np.floor((k+1)*mcx))])
        itschunks.append(parlist[np.int(np.floor((multiproc-1)*mcx)):nmc])

        plist = []
        for k in range(multiproc):
            p = Process(target=_formproc, args=(itschunks[k], pqueue))
            plist.append(p)
            p.start()

        for p in plist:
            p.join()

        ylist = []
        for k in range(multiproc):
            cylist = pqueue.get()
            ylist.extend(cylist)
            if verbose:
                estxy = np.average(np.array(ylist), axis=0)[0]
                print('mc:{0}/{1}: estxy[0]={2}'.
                      format(len(ylist), nmc, estxy))

        return ylist, np.average(np.array(ylist), axis=0), expvpara

    else:
        ylist = []
        if verbose:
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
                ylist.append((solfunc(cpar)).flatten())
            return ylist, np.average(np.array(ylist), axis=0), expvpara


def run_pce_sim_separable(solfunc=None, uncdims=None, abscissae=None,
                          multiproc=0):
    """ pce simulation for all PCE dimensions being the same
    """
    # compute the sols
    if multiproc > 1:
        pqueue = Queue()

        def comppart(itspart, partnum, queue):
            locylist = []
            for absctpl in itspart:
                locylist.append((solfunc(absctpl)).flatten())
            queue.put((partnum, locylist))

        lenits = abscissae.size**uncdims
        itspart = lenits/multiproc
        print('itspart {0} : lenits {1}'.format(itspart, lenits))
        itschunks = []
        for k in range(multiproc-1):
            itschunks.append(islice(product(abscissae, repeat=uncdims),
                             np.int(np.floor(k*itspart)),
                             np.int(np.floor((k+1)*itspart))))
        itschunks.append(islice(product(abscissae, repeat=uncdims),
                         np.int(np.floor((multiproc-1)*itspart)), lenits))
        plist = []
        for k in range(multiproc):
            p = Process(target=comppart, args=(itschunks[k], k, pqueue))
            plist.append(p)
            p.start()

        for p in plist:
            p.join()

        thinglist, orderarray = [], np.zeros((multiproc, ))
        for kk in range(multiproc):
            thething = pqueue.get()
            thinglist.append(thething[1])
            orderarray[thething[0]] = kk

        pceylist = []
        for kk in orderarray:
            pceylist.extend(thinglist[np.int(kk)])

    else:
        pceylist = []
        for absctpl in product(abscissae, repeat=uncdims):
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
        for idx, wtpl in enumerate(product(weights, repeat=uncdims)):
            cw = (np.array(wtpl)).prod()
            expv += cw*maty[:, idx]
        return scalefac*expv

    return abscissae, weights, comp_expv
