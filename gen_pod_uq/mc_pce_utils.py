from itertools import product, islice
from multiprocessing import Process
import os

import numpy as np
from scipy.io import savemat, loadmat

import multidim_galerkin_pod.chaos_expansion_utils as ceu


def solfunc_to_variance(solfunc, expv):
    def sfvariance(para):
        return (solfunc(para) - expv)**2
    return sfvariance


def run_mc_sim(parlist, solfunc, chunks=12, multiproc=0,
               tmpddir='_tmp_data_chunks/',
               comp_para_ev=True, verbose=False, ret_ev=True):
    expvpara = None
    if comp_para_ev:
        expvpara = np.average(np.array(parlist), axis=0)

    def _compallsols(paralist):
        locylist = []
        for cpara in paralist:
            locylist.append(solfunc(cpara).flatten())
        return locylist

    def _formproc(paralist, filestr):
        ylist = _compallsols(paralist)
        savemat(filestr, dict(ypart=np.array(ylist)))

    nmc = len(parlist)
    if multiproc > 1:
        mpid = os.getpid()
        mcx = nmc/multiproc

        itschunks = []
        fstrl = [tmpddir + f'_tmp_mc_chunk{1}of{multiproc}pid{mpid}']
        for k in range(multiproc-1):
            filestr = tmpddir + f'_tmp_mc_chunk{k+2}of{multiproc}pid{mpid}'
            fstrl.append(filestr)
            itschunks.append(parlist[np.int(np.floor(k*mcx)):
                                     np.int(np.floor((k+1)*mcx))])
        itschunks.append(parlist[np.int(np.floor((multiproc-1)*mcx)):nmc])

        plist = []
        for k, fstr in enumerate(fstrl):
            p = Process(target=_formproc, args=(itschunks[k], fstr))
            plist.append(p)
            p.start()

        for p in plist:
            p.join()

        yarray = loadmat(fstrl[0])['ypart']
        for fstr in fstrl[1:]:
            cychunk = loadmat(fstr)['ypart']
            yarray = np.vstack([yarray, cychunk])
            if verbose:
                estxy = np.average(yarray, axis=0)[0]
                print('mc:{0}/{1}: estxy[0]={2}'.
                      format(yarray.shape[0], nmc, estxy))

        return yarray, np.average(yarray, axis=0), expvpara

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

            return (np.array(ylist), np.average(np.array(ylist), axis=0),
                    expvpara)

        else:
            for cpar in parlist:
                ylist.append((solfunc(cpar)).flatten())
            return (np.array(ylist), np.average(np.array(ylist), axis=0),
                    expvpara)


def run_pce_sim_separable(solfunc=None, uncdims=None, abscissae=None,
                          tmpddir='_tmp_data_chunks/', multiproc=0):
    """ pce simulation for all PCE dimensions being the same
    """
    # compute the sols
    if multiproc > 1:
        # pqueue = Queue()

        def comppart(itspart, filestr):
            locylist = []
            for absctpl in itspart:
                locylist.append((solfunc(absctpl)).flatten().tolist())
            savemat(filestr, dict(ypart=np.array(locylist)))

        lenits = abscissae.size**uncdims
        itspart = lenits/multiproc
        # print('itspart {0} : lenits {1}'.format(itspart, lenits))
        itschunks = []
        for k in range(multiproc-1):
            itschunks.append(islice(product(abscissae, repeat=uncdims),
                             np.int(np.floor(k*itspart)),
                             np.int(np.floor((k+1)*itspart))))
        itschunks.append(islice(product(abscissae, repeat=uncdims),
                         np.int(np.floor((multiproc-1)*itspart)), lenits))
        plist = []
        fstrl = []
        for k in range(multiproc):
            filestr = tmpddir + '_tmp_pce_chunk{0}of{1}'.format(k+1, multiproc)
            p = Process(target=comppart, args=(itschunks[k], filestr))
            plist.append(p)
            fstrl.append(filestr)
            p.start()

        for p in plist:
            p.join()
        # print(pqueue.get())
        # import ipdb
        # ipdb.set_trace()

        ychunkl = []
        for fstr in fstrl:
            cychunk = loadmat(fstr)['ypart']
            ychunkl.append(cychunk)

        yarray = np.vstack(ychunkl)
        ypcedims = [yarray.shape[1]]
        ypcedims.extend([len(abscissae)]*uncdims)
        ypcedims = tuple(ypcedims)
        # arrange it in the tensor
        # first dimension is the state, then comes the uncertainty
        yrslttns = (yarray.T).reshape(ypcedims)

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

    def comp_vrnc(ytens, expv):
        ydim = ytens.shape[0]
        vrnc = 0
        matysqrd = np.square(ytens.reshape((ydim, -1)))
    # for kk, cw in enumerate(weights):
        for idx, wtpl in enumerate(product(weights, repeat=uncdims)):
            cw = (np.array(wtpl)).prod()
            vrnc += cw*matysqrd[:, idx]
        return scalefac*vrnc - expv**2

    return abscissae, weights, comp_expv, comp_vrnc


def pce_comp_vrnc(ytens, expv, weights=None, uncdims=None, scalefac=None):
    if scalefac is None:
        scalefac = (1./weights.sum())**uncdims
    ydim = ytens.shape[0]
    vrnc = 0
    matysqrd = np.square(ytens.reshape((ydim, -1)))
    for idx, wtpl in enumerate(product(weights, repeat=uncdims)):
        cw = (np.array(wtpl)).prod()
        vrnc += cw*matysqrd[:, idx]
    return scalefac*vrnc - expv**2
