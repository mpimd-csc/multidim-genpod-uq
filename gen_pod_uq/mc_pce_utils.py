from itertools import product, islice

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
            ylist.append((solfunc(cpar)).flatten())
        return ylist, np.average(np.array(ylist), axis=0), expvpara


def run_pce_sim_separable(solfunc=None, uncdims=None, abscissae=None,
                          multiproc=4):
    """ pce simulation for all PCE dimensions being the same
    """
    # compute the sols
    if multiproc > 1:
        from multiprocessing import Process, Queue
        pqueue = Queue()

        def comppart(itspart, partnum, queue):
            locylist = []
            for absctpl in itspart:
                locylist.append((solfunc(absctpl)).flatten())
            queue.put((partnum, locylist))

        if multiproc == 2:
            print('we use 2 processes')
            lenits = abscissae.size**uncdims
            itspart = np.int(np.floor(lenits/2))
            partone = islice(product(abscissae, repeat=uncdims), 0, itspart)
            parttwo = islice(product(abscissae, repeat=uncdims),
                             itspart, lenits)

            pone = Process(target=comppart, args=(partone, 1, pqueue))
            ptwo = Process(target=comppart, args=(parttwo, 2, pqueue))
            pone.start()
            ptwo.start()
            pone.join()
            ptwo.join()
            thingone = pqueue.get()
            thingtwo = pqueue.get()

            if thingone[0] == 1:
                pceylist = thingone[1]
                pceylist.extend(thingtwo[1])
            else:
                pceylist = thingtwo[1]
                pceylist.extend(thingone[1])

        elif multiproc == 4:
            print('we use 4 processes')
            lenits = abscissae.size**uncdims
            itspart = np.int(np.floor(lenits/4))
            partone = islice(product(abscissae, repeat=uncdims), 0, itspart)
            parttwo = islice(product(abscissae, repeat=uncdims),
                             itspart, 2*itspart)
            partthr = islice(product(abscissae, repeat=uncdims),
                             2*itspart, 3*itspart)
            partfou = islice(product(abscissae, repeat=uncdims),
                             3*itspart, lenits)

            pone = Process(target=comppart, args=(partone, 1, pqueue))
            ptwo = Process(target=comppart, args=(parttwo, 2, pqueue))
            pthr = Process(target=comppart, args=(partthr, 3, pqueue))
            pfou = Process(target=comppart, args=(partfou, 4, pqueue))
            pone.start()
            ptwo.start()
            pthr.start()
            pfou.start()
            pone.join()
            ptwo.join()
            pthr.join()
            pfou.join()
            thinglist, orderarray = [], np.zeros((4, ))
            for kk in range(4):
                thething = pqueue.get()
                thinglist.append(thething[1])
                orderarray[thething[0]-1] = kk

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
