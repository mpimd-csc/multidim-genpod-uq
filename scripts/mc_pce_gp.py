import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt
import time
import copy
import json

import logging
from rich.logging import RichHandler

import dolfin

# import spacetime_galerkin_pod.chaos_expansion_utils as ceu
import multidim_galerkin_pod.ten_sor_utils as tsu
import multidim_galerkin_pod.gen_pod_utils as gpu
from multidim_galerkin_pod.ldfnp_ext_cholmod import SparseFactorMassmat

import gen_pod_uq.mc_pce_utils as mpu

import dolfin_navier_scipy.data_output_utils as dou

from circle_subsec import get_problem
from cyl_subsec import get_problem as cylinder


def simit(problem='circle', meshlevel=None,
          mcruns=None, pcedimlist=None, plotplease=False,
          mcplease=False, pceplease=False, rbplease=False,
          mcpod=False, pcepod=False,
          rbparams={},
          checkredmod=False, pcexpy=None, pcvrnc=None,
          mcxpy=None, redmcruns=None,
          mcsnap=None, onlymeshtest=False,
          trainpcedim=None, targetpcedim=None,
          # basenu=5e-4, varia=-1e-4, varib=1e-4,
          plotpcepoddiff=False, pcepoddiffdim=9,
          multiproc=0, timings=1,
          nulb=6e-4, nuub=8e-4,
          distribution='uniform',
          databasemoments='cached-data/computed-moments.json',
          basisfrom='pce', poddimlist=[5, 8, 12]):

    if problem == 'cylinder':
        (get_sol, get_output, problemfems, plotit,
         get_red_problem) = cylinder(meshlevel=meshlevel)
        uncdims = 4
    else:
        get_sol, get_output, problemfems, get_red_problem = get_problem()
        uncdims = 5

    getsample = mpu.get_nu_sample(distribution=distribution, uncdims=uncdims,
                                  nulb=nulb, nuub=nuub)

    logging.info('Problem dimension: {0}'.format(problemfems['mmat'].shape[0]))

    nua, nub = nulb, nuub
    basenu = .5*(nua+nub)
    logging.info(f'nu in [{nua}, {nub}], median: {basenu}')
    logging.info('Distributed as: ' + distribution)

    cmat = problemfems['cmat']
    mmat = problemfems['mmat']

    prbsetupstr = 'N{0}nu{1:.2e}--{2:.2e}'.format(meshlevel, nulb, nuub) \
        + distribution

    if pcepod:
        filestr = prbsetupstr + '_pcepod'
    elif mcpod:
        filestr = prbsetupstr + '_mcpod'
    else:
        filestr = prbsetupstr
    if basisfrom == 'pce':
        bssstr = basisfrom + '{0}'.format(trainpcedim)
    elif basisfrom == 'mc':
        bssstr = basisfrom + '{0}_runs{1}'.format(mcsnap, timings)
    elif basisfrom == 'rb':
        bssstr = basisfrom + '_' + rbparams['samplemethod'] + \
            '{0}_runs{1}'.format(rbparams['nsample'], timings)
    else:
        pass
    filestr = filestr + '_bf' + bssstr + '.json'

    if onlymeshtest or plotplease:
        basenulist = [basenu]*uncdims
        basev = get_sol(basenulist)
        logging.info('N{1}: y(basenu)={0}'.format(cmat.dot(basev), meshlevel))
        basepvdfile = dolfin.File('results/basesol-nu{1:0.2e}-N{0}.pvd'.
                                  format(meshlevel, basenu))
        plotit(vvec=basev, pvdfile=basepvdfile, plotplease=plotplease)
        if onlymeshtest:
            return problemfems['mmat'].shape[0], cmat.dot(basev)

    if rbplease or basisfrom == 'rb':

        def get_rbbas(nsamples=None, nrbvecs=None):
            if rbparams['samplemethod'] == 'random':
                # ## XXX: here is randomness
                rbtrainnu = getsample(nsamples)
                logging.info(f'RB with {nsamples} random training points')
            else:
                raise NotImplementedError()
            # ## BLUNT implementation of RB
            # We precompute all solutions at the training parameters
            # to find *the max* without an estimator
            # (at least, we can use this solutions than to setup the RB)
            logging.info('computing all values for the RB train set ...')
            rbtrainset, _, _ = mpu.run_mc_sim(rbtrainnu, get_sol, verbose=True,
                                              multiproc=multiproc)
            logging.info('... done!')

            def _getmaxparam(cp_get_sol, dffun):
                mxdiff, mxdfpara, mxdfprid = 0, None, None
                diffl = []
                for cprid, cpara in enumerate(rbtrainnu):
                    cdiff = dffun(cp_get_sol(cpara),
                                  rbtrainset[cprid].reshape((-1, 1)))
                    # print(f'err: {cdiff} -- para val {cpara}')
                    if cdiff > mxdiff:
                        logging.debug(f'n-max: err: {cdiff.flatten()[0]:2e}' +
                                      f' \n-- pv: {cpara}')
                        mxdfpara, mxdfprid = cpara, cprid
                        mxdiff = cdiff
                    diffl.append(cdiff)
                logging.debug(f'found max err: {mxdiff.flatten()[0]:2e} \n' +
                              f'-- at: {mxdfpara}')
                return mxdfpara, rbtrainset[mxdfprid].reshape((-1, 1))

            def _dffun(vone, vtwo):
                diffv = vone-vtwo
                return np.sqrt(diffv.T @ mmat @ diffv)

            rbbas = rbtrainset[0].reshape((-1, 1))

            for rbdim in range(nrbvecs-1):
                logging.debug('*** RB: Greedy iteration ' +
                              f'{rbdim+1}/{nrbvecs-1} ***')
                crb_red_realize_sol, _, _, _ = get_red_problem(rbbas)

                def _compun(para):
                    return rbbas @ crb_red_realize_sol(para)
                mxdfpara, mxdfsol = _getmaxparam(_compun, dffun=_dffun)
                # print('returned: ', mxdfpara)
                # rbcheck = get_sol(mxdfpara)
                # print('right vec?: ', dffun(rbcheck, mxdfsol))
                rbbas = np.hstack([rbbas, mxdfsol])

            return rbbas

        if rbplease:
            rbbas = get_rbbas(nsamples=rbparams['nsample'],
                              nrbvecs=rbparams['N'])

    # ## CHAP Monte Carlo
    if mcplease:
        # varinu = nulb + (nulb-varia)*np.random.rand(mcruns, uncdims)
        # ## XXX: here is randomness
        varinu = getsample(mcruns)
        expvnu = np.average(varinu, axis=0)
        logging.info('expected value of nu: ', expvnu)
        varinulist = varinu.tolist()
        mcout, mcxpy, expvnu = mpu.run_mc_sim(varinulist, get_output,
                                              verbose=True,
                                              multiproc=multiproc)

        mmcsolfile = dolfin.File('results/mmcsol.pvd')
        curv = get_sol(expvnu.tolist())
        plotit(vvec=curv, pvdfile=mmcsolfile, plotplease=plotplease)
        logging.info('y(estxnu)={0}'.format(cmat.dot(curv)))

        if plotplease:
            plt.figure(89)
            plt.plot(mcout, '.')
            plt.show()

    def put_mmnts_db(pcedim=None, expv=None, vrnc=None):

        with open(databasemoments, 'r+') as fjs:
            try:
                dbmmnts = json.load(fjs)
            except json.JSONDecodeError:
                dbmmnts = {}

        nukey = f'{nua:.1e}to{nub:.1e}'
        try:
            subdict = dbmmnts[f'{meshlevel}'][nukey][distribution]
        except KeyError:
            dbmmnts.update({f'{meshlevel}': {nukey: {distribution: {}}}})
            subdict = dbmmnts[f'{meshlevel}'][nukey][distribution]
        subdict.update({f'{pcedim}': {'expv': expv.flatten()[0],
                                      'vrnc': vrnc.flatten()[0]}})

        with open(databasemoments, 'w') as fjs:
            fjs.write(json.dumps(dbmmnts))

        return

    def get_mmnts_db(pcedim):
        with open(databasemoments, 'r+') as fjs:
            try:
                dbmmnts = json.load(fjs)
            except json.JSONDecodeError:
                dbmmnts = {}
        nukey = f'{nua:.1e}to{nub:.1e}'
        try:
            subdict = dbmmnts[f'{meshlevel}'][nukey][distribution]
            pcexpy = subdict[f'{pcedim}']['expv']
            pcvrnc = subdict[f'{pcedim}']['vrnc']
            return pcexpy, pcvrnc, True
        except KeyError as e:
            logging.debug(e, exc_info=True)
            return None, None, False

    # ## CHAP Polynomial Chaos Expansion
    if pceplease:
        for pcedim in pcedimlist:
            pcexpy, pcvrnc, esth = get_mmnts_db(pcedim=pcedim)
            if esth:
                logging.info('PCE({0}): E(y): {1}'.format(pcedim, pcexpy))
                logging.info('PCE({0}): V(y): {1}'.format(pcedim, pcvrnc))
                logging.debug('loaded from ' + databasemoments)
            else:
                logging.info(f'Computing: PCE({pcedim})')

                # ## XXX: here is randomness
                abscissae, weights, compexpv, _ = mpu.\
                    setup_pce(distribution=distribution,
                              distrpars=dict(a=nua, b=nub),
                              pcedim=pcedim, uncdims=uncdims)
                ysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
                                                     uncdims=uncdims,
                                                     multiproc=multiproc,
                                                     abscissae=abscissae)
                pcexpy = compexpv(ysoltens)
                pcexpysqrd = compexpv(np.square(ysoltens))
                pcvrnc = pcexpysqrd-pcexpy**2
                logging.info('PCE({0}): E(y): {1}'.format(pcedim, pcexpy))
                logging.info('PCE({0}): E(yy): {1}'.format(pcedim, pcexpysqrd))
                logging.info('PCE({0}): V(y): {1}'.format(pcedim, pcvrnc))
                put_mmnts_db(pcedim=pcedim, expv=pcexpy, vrnc=pcvrnc)

    # if rbplease:
    #     rbey = np.mean(cmat @ rbbas)
    #     print(f'RB({rbparams["N"]}): E(y): {rbey}')

    if plotpcepoddiff:
        pcedim = pcedimlist[-1]
        pcepoddiffstr = 'pcepoddiff{0}_'.format(pcedim) + filestr
        try:
            pxexpxdct = dou.load_json_dicts(pcepoddiffstr)
            pcexpx = np.array(pxexpxdct['pcexpx'])
            logging.info('loaded the pce-Ex from: ', pcepoddiffstr)
        except IOError:
            # ## XXX: here is randomness
            abscissae, weights, compexpv, _ = mpu.\
                setup_pce(distribution=distribution,
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcedim, uncdims=uncdims)
            xsoltens = mpu.run_pce_sim_separable(solfunc=get_sol,
                                                 uncdims=uncdims,
                                                 multiproc=multiproc,
                                                 abscissae=abscissae)
            pcexpx = compexpv(xsoltens)
            jsfile = open(pcepoddiffstr, mode='w')
            jsfile.write(json.dumps({'pcexpx': pcexpx.tolist(),
                                     'podpcexpx': {}}))
            jsfile.close()
            logging.info('saved the pce-Ex to: ', pcepoddiffstr)

    if not (pcepod or mcpod):
        return

    # ## CHAP genpod
    mmat = problemfems['mmat']
    facmy = SparseFactorMassmat(mmat)

    truthexpy, truthvrnc, esth = get_mmnts_db(pcedim=pcedimlist[-1])
    if esth:
        logging.info(f'truth: PCE({pcedimlist[-1]}): E(y): {truthexpy}')
        logging.info(f'truth: PCE({pcedimlist[-1]}): V(y): {truthvrnc}')
        logging.debug('loaded from ' + databasemoments)
    else:
        logging.info('no truthvalues loaded -- tried with ' + databasemoments)
        logging.debug(f'for with pcedim={pcedimlist[-1]}')
    tdict = dict(truthvrnc=truthvrnc, truthexpy=truthexpy)

    np.random.seed(1)  # seed for the random `mc` basis

    for tit in range(timings):
        loctdict = {'basisfrom': basisfrom}
        if basisfrom == 'pce':
            trttstart = time.time()
            # ## XXX: here is randomness
            trnabscissae, trnweights, trncompexpv, trncomvrnc = mpu.\
                setup_pce(distribution=distribution,
                          distrpars=dict(a=nua, b=nub),
                          pcedim=trainpcedim, uncdims=uncdims)
            pcewmatfac = sps.dia_matrix((np.sqrt(trnweights), 0),
                                        shape=(trainpcedim, trainpcedim))

            # the multidim mass matrices
            mfl = [facmy.F]
            mfl.extend([pcewmatfac]*uncdims)

            logging.info('Snapshot computation started...')
            trainsoltens = mpu.run_pce_sim_separable(solfunc=get_sol,
                                                     uncdims=uncdims,
                                                     multiproc=multiproc,
                                                     abscissae=trnabscissae)
            # cysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
            #                                       uncdims=uncdims,
            #                                       abscissae=abscissae)
            trtelt = time.time() - trttstart
            logging.info(f'Snapshot computation: Elapsed time: {trtelt}')
            trainexpv = trncompexpv(trainsoltens)
            trainpcexpy = cmat.dot(trainexpv)
            logging.info(f'estimated expected value (pce): {trainpcexpy}')
            loctdict.update({'training-pce-expv': trainpcexpy.tolist(),
                             'traintime': trtelt})

            if pcexpy is not None:
                trnrpcexpy = (trainpcexpy-pcexpy)
                logging.info(
                    '-> difference expv (pce): {0}'.format(trnrpcexpy))
            if mcxpy is not None:
                trnrmcexpy = mcxpy - trainpcexpy
                logging.info(
                    '-> difference mc estimate: {0}'.format(trnrmcexpy))

            def get_pod_vecs(poddim=None):
                return tsu.modeone_massmats_svd(trainsoltens, mfl, poddim)

        elif basisfrom == 'mc':
            trttstart = time.time()
            # ## XXX: here is randomness
            varinu = getsample(mcsnap)
            expvnu = np.average(varinu, axis=0)
            varinulist = varinu.tolist()
            mcout, _, _ = mpu.run_mc_sim(varinulist, get_sol,
                                         multiproc=multiproc)
            lymcmat = facmy.Ft*mcout.T
            trtelt = time.time() - trttstart
            logging.info('POD basis by {0} random samplings'.format(mcsnap))
            snpshmean = np.average(mcout.T, axis=1)
            snpshymean = cmat.dot(snpshmean)
            logging.info('estimated mean of the samplings: {snpshymean}')
            loctdict.update({'training-mc-estmean': snpshymean.tolist(),
                             'traintime': trtelt})
            if pcexpy is not None:
                trnrpcexpy = pcexpy - snpshmean
                logging.info(
                    '-> difference expv (pce): {0}'.format(trnrpcexpy))
            if mcxpy is not None:
                trnrmcexpy = mcxpy - np.average(cmat.dot(mcout.T), axis=1)
                logging.info(
                    '-> difference mc estimate: {0}'.format(trnrmcexpy))

            def get_pod_vecs(poddim=None):
                ypodvecs = gpu.get_ksvvecs(sol=lymcmat, poddim=poddim,
                                           plotsvs=plotplease, labl='SVs')
                return ypodvecs

        elif basisfrom == 'rb':
            trttstart = time.time()
            rbbas = get_rbbas(nsamples=rbparams['nsample'],
                              nrbvecs=poddimlist[-1])
            trtelt = time.time() - trttstart
            logging.info(f"RB basis from {rbparams['nsample']} random samples")
            lymrbvecs = facmy.Ft*rbbas

            loctdict.update({'traintime': trtelt})

            def get_pod_vecs(poddim=None):
                ''' return the first `poddim` RB vectors

                orthogonal wrt `M` inner product

                to have the reduction and projection defined
                in line with `pcepod`
                '''
                ypodvecs = gpu.get_ksvvecs(sol=lymrbvecs[:, :poddim],
                                           poddim=poddim,
                                           plotsvs=plotplease, labl='SVs')
                # no reduction, just orthogonalization
                return ypodvecs

        else:
            raise NotImplementedError()

        # lypceymat = pceymat
        redsolfile = dolfin.File('results/rdsol-N{0}pod.pvd'.format(meshlevel))

        pcepoddict = {}
        mcpoddict = {}
        crmeltlist = []
        rmprjerrs = []
        for poddim in poddimlist:
            tstart = time.time()
            ypodvecs = get_pod_vecs(poddim)
            # if basisfrom == 'rb':
            #     lyitVy = ypodvecs
            # else:
            lyitVy = facmy.solve_Ft(ypodvecs)
            red_realize_sol, red_realize_output, red_probfems, red_plotit \
                = get_red_problem(lyitVy)
            red_cmat = red_probfems['cmat']
            crmelt = time.time() - tstart
            crmeltlist.append(crmelt)

            logging.info(f'poddim:{poddim}: ctime ROM: elt: {crmelt}')
            if basisfrom == 'pce':  # or basisfrom == 'rb':
                cndsdexpv = lyitVy.T.dot(mmat.dot(trainexpv))
                prjerror = trainpcexpy - red_cmat.dot(cndsdexpv)
                rmprjerrs.append(prjerror.tolist())
            elif basisfrom == 'mc':
                cndsshm = lyitVy.T.dot(mmat.dot(snpshmean))
                prjerror = snpshymean - red_cmat.dot(cndsshm)
                rmprjerrs.append(prjerror.tolist())
            else:
                pass

            if checkredmod:
                nulist = [basenu]*uncdims
                redv = red_realize_sol(nulist)
                red_plotit(vvec=redv, pvdfile=redsolfile,
                           plotplease=plotplease)
                redy = red_cmat.dot(redv)
                logging.info(f'N{meshlevel}pod{poddim}red_y(basenu)={redy}')

            if pcepod:
                pcereslist, pcepodeysqrd, eltlist = [], [], []
                logging.info('dim of reduced model: {0}'.format(poddim))
                for pcedim in pcedimlist:
                    # ## XXX: here is randomness
                    abscissae, weights, compredexpv, compredvrnc = mpu.\
                        setup_pce(distribution=distribution,
                                  distrpars=dict(a=nua, b=nub),
                                  pcedim=pcedim, uncdims=uncdims)
                    tstart = time.time()
                    redysoltens = mpu.\
                        run_pce_sim_separable(solfunc=red_realize_output,
                                              multiproc=multiproc,
                                              uncdims=uncdims,
                                              abscissae=abscissae)
                    redpcexpy = compredexpv(redysoltens)
                    elt = time.time() - tstart
                    redpcexpeysqrd = compredexpv(np.square(redysoltens))
                    pcereslist.append(redpcexpy.tolist())
                    pcepodeysqrd.append(redpcexpeysqrd.tolist())
                    eltlist.append(elt)
                    if truthexpy is not None:
                        logging.info(f'pce={pcedim:2.0f}: ' + f'elt={elt:.2f}')
                        logging.info(f'e_xpvl={(redpcexpy-truthexpy)[0]:.3e}')
                        if truthvrnc is not None:
                            evrnc = redpcexpeysqrd - redpcexpy**2 - truthvrnc
                            logging.info(f'e_vrnc={evrnc[0]:.3e}')

                pcepoddict.update({poddim: {'pcedims': pcedimlist,
                                            'pceres': pcereslist,
                                            'pcepodeyys': pcepodeysqrd,
                                            'elts': eltlist}})
            if mcpod:
                varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
                expvnu = np.average(varinu, axis=0)
                logging.info('expected value of nu: ', expvnu)
                varinulist = varinu.tolist()
                mcptstart = time.time()
                (mcout, rmcxpy,
                 expvnu) = mpu.run_mc_sim(varinulist, red_realize_output,
                                          multiproc=multiproc)
                mcpelt = time.time() - mcptstart
                if mcxpy is not None:
                    ifs = f'mcruns={0:2.0f}, poddim={poddim:2.0f}, ' + \
                        'rmcxpy-mcxpy={0}'.format(rmcxpy-mcxpy)
                    logging.info(ifs)
                mcpoddict.update({poddim: {'mcruns': mcruns,
                                           'mcres': rmcxpy.tolist(),
                                           'elt': mcpelt}})
        if pcepod:
            loctdict.update({'pcepod': copy.deepcopy(pcepoddict)})
        if mcpod:
            loctdict.update({'mcpod': copy.deepcopy(mcpoddict)})
        loctdict.update({'comp-redmod-elts': crmeltlist,
                         'redmod-prj-errs': rmprjerrs})

        tdict.update({tit: copy.deepcopy(loctdict)})

    jsfile = open(filestr, mode='w')
    jsfile.write(json.dumps(tdict))
    logging.info('output saved to ' + filestr)

    if plotpcepoddiff:
        pxexpxdct = dou.load_json_dicts(pcepoddiffstr)
        pcexpx = np.array(pxexpxdct['pcexpx'])
        try:
            podpcexpx = np.array(pxexpxdct['podpcexpx'][pcepoddiffdim])
            logging.info(f'loaded the pod{pcepoddiffdim}-pce-Ex from: ' +
                         pcepoddiffstr)
        except KeyError:
            ypodvecs = get_pod_vecs(pcepoddiffdim)
            lyitVy = facmy.solve_Ft(ypodvecs)
            red_realize_sol, red_realize_output, red_probfems, red_plotit \
                = get_red_problem(lyitVy)
            red_cmat = red_probfems['cmat']
            # ## XXX: here is randomness
            abscissae, weights, compredexpv, compredvrnc = mpu.\
                setup_pce(distribution=distribution,
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcedimlist[-1], uncdims=uncdims)
            redxsoltens = mpu.\
                run_pce_sim_separable(solfunc=red_realize_sol,
                                      multiproc=multiproc,
                                      uncdims=uncdims,
                                      abscissae=abscissae)
            podpcexpx = compredexpv(redxsoltens)
            pxexpxdct['podpcexpx'].update({pcepoddiffdim: podpcexpx.tolist()})
            logging.info(f'appended the pod{pcepoddiffdim}-pce-Ex to: ' +
                         pcepoddiffstr)
        ppdpvdfile = dolfin.File('results/pce{2}pod{3}dif-nu{1:0.2e}-N{0}.pvd'.
                                 format(meshlevel, basenu, pcedimlist[-1],
                                        pcepoddiffdim))
        plotit(vvec=np.atleast_2d(pcexpx-lyitVy.dot(podpcexpx)).T,
               pvdfile=ppdpvdfile, plotplease=True)

    plt.show()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, handlers=[RichHandler()],
                        format='%(message)s',
                        datefmt="[%X]",
                        )
    problem = 'cylinder'
    distribution = 'beta-2-5'
    distribution = 'uniform'
    meshlevel = 4
    mcruns = 10  # 200
    pcedimlist = [2, 4, 5]  # , 3, 4]  # , 3, 4, 5]  # , 7]
    multiproc = 4
    timings = 1
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
    basisfrom = 'rb'
    basisfrom = 'pce'
    rbparams = dict(samplemethod='random', nsample=16, N=16)

    simit(mcruns=mcruns, pcedimlist=pcedimlist, problem=problem,
          distribution=distribution,
          meshlevel=meshlevel, timings=timings,
          plotplease=plotplease, basisfrom=basisfrom, multiproc=multiproc,
          rbparams=rbparams, trainpcedim=2, targetpcedim=5,
          mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
