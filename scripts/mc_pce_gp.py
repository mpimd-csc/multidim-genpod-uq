import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt
import time
import copy
import json

import dolfin

# import spacetime_galerkin_pod.chaos_expansion_utils as ceu
import multidim_galerkin_pod.ten_sor_utils as tsu
import multidim_galerkin_pod.gen_pod_utils as gpu
from multidim_galerkin_pod.ldfnp_ext_cholmod import SparseFactorMassmat

import gen_pod_uq.mc_pce_utils as mpu

import dolfin_navier_scipy.data_output_utils as dou

from circle_subsec import get_problem
from cyl_subsec import get_problem as cylinder

plotpcepoddiff = True
pcepoddiffdim = 9


def simit(problem='circle', meshlevel=None,
          mcruns=None, pcedimlist=None, plotplease=False,
          mcplease=False, pceplease=False, mcpod=False, pcepod=False,
          checkredmod=False, pcexpy=None, pcvrnc=0.,
          mcxpy=None, redmcruns=None,
          mcsnap=None, pcesnapdim=None, onlymeshtest=False,
          # basenu=5e-4, varia=-1e-4, varib=1e-4,
          multiproc=0, timings=1,
          nulb=6e-4, nuub=8e-4,
          basisfrom='pce', poddimlist=[5, 10, 20]):

    if problem == 'cylinder':
        (get_sol, get_output, problemfems, plotit,
         get_red_problem) = cylinder(meshlevel=meshlevel)
        uncdims = 4
    else:
        get_sol, get_output, problemfems, get_red_problem = get_problem()
        uncdims = 5

    print(problemfems['mmat'].shape[0])

    nua, nub = nulb, nuub
    basenu = .5*(nua+nub)
    print('basenu: {0}'.format(basenu))
    print('mininu: {0}'.format(nua))
    print('maxinu: {0}'.format(nub))

    cmat = problemfems['cmat']

    filestr = 'N{0}nu{1:.2e}--{2:.2e}'.format(meshlevel, nulb, nuub)
    if pcepod:
        filestr = filestr + '_pcepod{0}'.format(pcesnapdim)
    if mcpod:
        filestr = filestr + '_mcpod{0}'.format(mcsnap)
    filestr = filestr + '_bf' + basisfrom + '.json'

    if onlymeshtest or plotplease:
        basenulist = [basenu]*uncdims
        basev = get_sol(basenulist)
        print('N{1}: y(basenu)={0}'.format(cmat.dot(basev), meshlevel))
        basepvdfile = dolfin.File('results/basesol-nu{1:0.2e}-N{0}.pvd'.
                                  format(meshlevel, basenu))
        plotit(vvec=basev, pvdfile=basepvdfile, plotplease=plotplease)
        if onlymeshtest:
            return problemfems['mmat'].shape[0], cmat.dot(basev)

    # ## CHAP Monte Carlo
    if mcplease:
        # varinu = nulb + (nulb-varia)*np.random.rand(mcruns, uncdims)
        varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
        expvnu = np.average(varinu, axis=0)
        print('expected value of nu: ', expvnu)
        varinulist = varinu.tolist()
        mcout, mcxpy, expvnu = mpu.run_mc_sim(varinulist, get_output,
                                              verbose=True,
                                              multiproc=multiproc)

        mmcsolfile = dolfin.File('results/mmcsol.pvd')
        curv = get_sol(expvnu.tolist())
        plotit(vvec=curv, pvdfile=mmcsolfile, plotplease=plotplease)
        print('y(estxnu)={0}'.format(cmat.dot(curv)))

        if plotplease:
            plt.figure(89)
            plt.plot(mcout, '.')
            plt.show()

    # ## CHAP Polynomial Chaos Expansion
    if pceplease:
        for pcedim in pcedimlist:
            abscissae, weights, compexpv, _ = mpu.\
                setup_pce(distribution='uniform',
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcedim, uncdims=uncdims)
            # abscarray, weightsarray = np.array(abscissae), np.array(weights)
            ysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
                                                 uncdims=uncdims,
                                                 multiproc=multiproc,
                                                 abscissae=abscissae)
            pcexpy = compexpv(ysoltens)
            pcexpysqrd = compexpv(np.square(ysoltens))
            print('PCE({0}): E(y): {1}'.format(pcedim, pcexpy))
            print('PCE({0}): V(y): {1}'.format(pcedim, pcexpysqrd-pcexpy**2))

    if plotpcepoddiff:
        pcedim = pcedimlist[-1]
        pcepoddiffstr = 'pcepoddiff{0}_'.format(pcedim) + filestr
        try:
            pxexpxdct = dou.load_json_dicts(pcepoddiffstr)
            pcexpx = np.array(pxexpxdct['pcexpx'])
        except IOError:
            abscissae, weights, compexpv, _ = mpu.\
                setup_pce(distribution='uniform',
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

    if not (pcepod or mcpod):
        return

    # ## CHAP genpod
    mmat = problemfems['mmat']
    facmy = SparseFactorMassmat(mmat)

    tdict = {}

    np.random.seed(1)  # seed for the random `mc` basis

    for tit in range(timings):
        loctdict = {'basisfrom': basisfrom}
        if basisfrom == 'pce':
            trttstart = time.time()
            trnabscissae, trnweights, trncompexpv, trncomvrnc = mpu.\
                setup_pce(distribution='uniform',
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcesnapdim, uncdims=uncdims)
            pcewmatfac = sps.dia_matrix((np.sqrt(trnweights), 0),
                                        shape=(pcesnapdim, pcesnapdim))

            mfl = [facmy.F]
            mfl.extend([pcewmatfac]*uncdims)
            trainsoltens = mpu.run_pce_sim_separable(solfunc=get_sol,
                                                     uncdims=uncdims,
                                                     multiproc=multiproc,
                                                     abscissae=trnabscissae)
            # cysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
            #                                       uncdims=uncdims,
            #                                       abscissae=abscissae)
            trtelt = time.time() - trttstart
            print('{0}: Elapsed time: {1}'.format('snapshot computation',
                                                  trtelt))
            trainexpv = trncompexpv(trainsoltens)
            trainpcexpy = cmat.dot(trainexpv)
            print('estimated expected value (pce): {0}'.format(trainpcexpy))
            loctdict.update({'training-pce-expv': trainpcexpy.tolist(),
                             'traintime': trtelt})

            if pcexpy is not None:
                trnrpcexpy = (trainpcexpy-pcexpy)
                print('-> difference expv (pce): {0}'.format(trnrpcexpy))
            if mcxpy is not None:
                trnrmcexpy = mcxpy - trainpcexpy
                print('-> difference mc estimate: {0}'.format(trnrmcexpy))

            def get_pod_vecs(poddim=None):
                return tsu.modeone_massmats_svd(trainsoltens, mfl, poddim)

        elif basisfrom == 'mc':
            trttstart = time.time()
            varinu = nulb + (nuub-nulb)*np.random.rand(mcsnap, uncdims)
            expvnu = np.average(varinu, axis=0)
            varinulist = varinu.tolist()
            mcout, _, _ = mpu.run_mc_sim(varinulist, get_sol,
                                         multiproc=multiproc)
            lymcmat = facmy.Ft*mcout.T
            trtelt = time.time() - trttstart
            print('POD basis by {0} random samplings'.format(mcsnap))
            snpshmean = np.average(mcout.T, axis=1)
            snpshymean = cmat.dot(snpshmean)
            print('estimated mean of the samplings: {0}'.format(snpshymean))
            loctdict.update({'training-mc-estmean': snpshymean.tolist(),
                             'traintime': trtelt})
            if pcexpy is not None:
                trnrpcexpy = pcexpy - snpshmean
                print('-> difference expv (pce): {0}'.format(trnrpcexpy))
            if mcxpy is not None:
                trnrmcexpy = mcxpy - np.average(cmat.dot(mcout.T), axis=1)
                print('-> difference mc estimate: {0}'.format(trnrmcexpy))

            def get_pod_vecs(poddim=None):
                ypodvecs = gpu.get_ksvvecs(sol=lymcmat, poddim=poddim,
                                           plotsvs=plotplease, labl='SVs')
                return ypodvecs

        # lypceymat = pceymat
        redsolfile = dolfin.File('results/rdsol-N{0}pod.pvd'.format(meshlevel))

        pcepoddict = {}
        mcpoddict = {}
        crmeltlist = []
        rmprjerrs = []
        for poddim in poddimlist:
            tstart = time.time()
            ypodvecs = get_pod_vecs(poddim)
            lyitVy = facmy.solve_Ft(ypodvecs)
            red_realize_sol, red_realize_output, red_probfems, red_plotit \
                = get_red_problem(lyitVy)
            red_cmat = red_probfems['cmat']
            crmelt = time.time() - tstart
            crmeltlist.append(crmelt)

            print('poddim:{2}: {0}: elt: {1}'.format('reduced model comp',
                                                     crmelt, poddim))
            if basisfrom == 'pce':
                cndsdexpv = lyitVy.T.dot(mmat.dot(trainexpv))
                prjerror = trainpcexpy - red_cmat.dot(cndsdexpv)
            elif basisfrom == 'mc':
                cndsshm = lyitVy.T.dot(mmat.dot(snpshmean))
                prjerror = snpshymean - red_cmat.dot(cndsshm)

            rmprjerrs.append(prjerror.tolist())

            if checkredmod:
                nulist = [basenu]*uncdims
                redv = red_realize_sol(nulist)
                red_plotit(vvec=redv, pvdfile=redsolfile,
                           plotplease=plotplease)
                print('N{1}pod{2}red_y(basenu)={0}'.format(red_cmat.dot(redv),
                                                           meshlevel, poddim))

            if pcepod:
                pcereslist, pcepodeysqrd, eltlist = [], [], []
                print('dim of reduced model: {0}'.format(poddim))
                for pcedim in pcedimlist:
                    abscissae, weights, compredexpv, compredvrnc = mpu.\
                        setup_pce(distribution='uniform',
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
                    if pcexpy is not None:
                        print('pce={0:2.0f}, exypce={1}, elt={2:.2f}'.
                              format(pcedim, redpcexpy-pcexpy, elt))
                    if pcvrnc is not None:
                        print('pce={0:2.0f}, evrnc={1}'.
                              format(pcedim, redpcexpeysqrd-pcexpy**2-pcvrnc))

                pcepoddict.update({poddim: {'pcedims': pcedimlist,
                                            'pceres': pcereslist,
                                            'pcepodeyys': pcepodeysqrd,
                                            'elts': eltlist}})

            if mcpod:
                varinu = nulb + (nuub-nulb)*np.random.rand(mcruns, uncdims)
                expvnu = np.average(varinu, axis=0)
                print('expected value of nu: ', expvnu)
                varinulist = varinu.tolist()
                mcptstart = time.time()
                (mcout, rmcxpy,
                 expvnu) = mpu.run_mc_sim(varinulist, red_realize_output,
                                          multiproc=multiproc)
                mcpelt = time.time() - mcptstart
                if mcxpy is not None:
                    print('mcruns={0:2.0f}, poddim={2:2.0f}, rmcxpy-mcxpy={1}'.
                          format(redmcruns, rmcxpy-mcxpy, poddim))
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
    print('output saved to ' + filestr)

    if plotpcepoddiff:
        pxexpxdct = dou.load_json_dicts(pcepoddiffstr)
        pcexpx = np.array(pxexpxdct['pcexpx'])
        try:
            podpcexpx = np.array(pxexpxdct['podpcexpx'][pcepoddiffdim])
        except KeyError:
            ypodvecs = get_pod_vecs(pcepoddiffdim)
            lyitVy = facmy.solve_Ft(ypodvecs)
            red_realize_sol, red_realize_output, red_probfems, red_plotit \
                = get_red_problem(lyitVy)
            red_cmat = red_probfems['cmat']
            abscissae, weights, compredexpv, compredvrnc = mpu.\
                setup_pce(distribution='uniform',
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcedimlist[-1], uncdims=uncdims)
            redxsoltens = mpu.\
                run_pce_sim_separable(solfunc=red_realize_sol,
                                      multiproc=multiproc,
                                      uncdims=uncdims,
                                      abscissae=abscissae)
            podpcexpx = compredexpv(redxsoltens)
            pxexpxdct['podpcexpx'].update({pcepoddiffdim: podpcexpx.tolist()})

    plt.show()


if __name__ == '__main__':
    problem = 'cylinder'
    meshlevel = 6
    mcruns = 10  # 200
    pcedimlist = [2, 3, 4]  # , 3, 4, 5]  # , 7]
    mcplease = False
    pceplease = False
    plotplease = False
    mcpod = False
    pcepod = False
    # ## make it come true
    # mcplease = True
    pceplease = True
    plotplease = True
    # pcepod = True
    # mcpod = True
    basisfrom = 'mc'
    basisfrom = 'pce'
    simit(mcruns=mcruns, pcedimlist=pcedimlist, problem=problem,
          meshlevel=meshlevel,
          plotplease=plotplease, basisfrom=basisfrom,
          mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
