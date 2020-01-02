import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

import dolfin

import spacetime_galerkin_pod.chaos_expansion_utils as ceu
import spacetime_galerkin_pod.ten_sor_utils as tsu
import spacetime_galerkin_pod.gen_pod_utils as gpu
from spacetime_galerkin_pod.ldfnp_ext_cholmod import SparseFactorMassmat

import gen_pod_uq.mc_pce_utils as mpu

from circle_subsec import get_problem
from cyl_subsec import get_problem as cylinder


def simit(problem='circle', meshlevel=None,
          mcruns=None, pcedimlist=None, plotplease=False,
          mcplease=False, pceplease=False, mcpod=False, pcepod=False,
          checkredmod=False, pcexpy=None, mcxpy=None, redmcruns=None,
          mcsnap=None, pcesnapdim=None,
          basisfrom='pce'):

    if problem == 'cylinder':
        (get_sol, get_output, problemfems, plotit,
         get_red_problem) = cylinder(meshlevel=meshlevel)
        uncdims = 4
    else:
        get_sol, get_output, problemfems, get_red_problem = get_problem()
        uncdims = 5

    print(problemfems['mmat'].shape[0])

    basenu = 1e-3
    varia = 0.
    varib = 5e-4
    nua, nub = basenu+varia, basenu+varib
    cmat = problemfems['cmat']

    basepvdfile = dolfin.File('results/basesol-N{0}.pvd'.format(meshlevel))
    basenulist = [basenu]*uncdims
    basev = get_sol(basenulist)
    plotit(vvec=basev, pvdfile=basepvdfile, plotplease=plotplease)

    print('N{1}: y(basenu)={0}'.format(cmat.dot(basev), meshlevel))

    # ## CHAP Monte Carlo
    if mcplease:
        varinu = basenu + (varib-varia)*np.random.rand(mcruns, uncdims)
        expvnu = np.average(varinu, axis=0)
        print('expected value of nu: ', expvnu)
        varinulist = varinu.tolist()
        mcout, mcxpy, expvnu = mpu.run_mc_sim(varinulist, get_output,
                                              verbose=True)

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
            abscissae, weights, compexpv = mpu.\
                setup_pce(distribution='uniform',
                          distrpars=dict(a=nua, b=nub),
                          pcedim=pcedim, uncdims=uncdims)
            # abscarray, weightsarray = np.array(abscissae), np.array(weights)
            ysoltens = mpu.run_pce_sim_separable(solfunc=get_output,
                                                 uncdims=uncdims,
                                                 abscissae=abscissae)
            pcexpy = compexpv(ysoltens)
            print('PCE({0}): E(y): {1}'.format(pcedim, pcexpy))

    if not (pcepod or mcpod):
        return

    # ## CHAP genpod
    pcedim = pcesnapdim
    mmat = problemfems['mmat']

    abscissae, weights = ceu.get_gaussqr_uniform(N=pcedim, a=nua, b=nub)
    facmy = SparseFactorMassmat(mmat)

    poddimlist = [5, 10, 20]  # , 40]

    # pcewmat = sps.dia_matrix((weights, 0), shape=(pcedim, pcedim))
    pcewmatfac = sps.dia_matrix((np.sqrt(weights), 0), shape=(pcedim, pcedim))

    mfl = [facmy.F]
    mfl.extend([pcewmatfac]*uncdims)

    if basisfrom == 'pce':
        ysoltens = mpu.run_pce_sim_separable(solfunc=get_sol,
                                             uncdims=uncdims,
                                             abscissae=abscissae)

        def get_pod_vecs(poddim=None):
            return tsu.modeone_massmats_svd(ysoltens, mfl, poddim)

    elif basisfrom == 'mc':
        varinu = basenu + (varib-varia)*np.random.rand(mcsnap, uncdims)
        expvnu = np.average(varinu, axis=0)
        varinulist = varinu.tolist()
        mcout, _, _ = mpu.run_mc_sim(varinulist, get_sol)
        pceymat = np.array(mcout).T
        lypceymat = facmy.Ft*pceymat
        print('POD basis by {0} random samplings'.format(mcsnap))

        def get_pod_vecs(poddim=None):
            ypodvecs = gpu.get_ksvvecs(sol=lypceymat, poddim=poddim,
                                       plotsvs=plotplease, labl='SVs')
            return ypodvecs

    # lypceymat = pceymat
    redsolfile = dolfin.File('results/redsol-N{0}pods.pvd'.format(meshlevel))
    rmcxpyl, rpcesxpyl = [], []
    cmpwmc = 'comp reduced mc and full mc'
    cmpwpce = 'comp reduced mc and full pce'

    for poddim in poddimlist:
        print('pcexpy: {0}'.format(pcexpy))
        ypodvecs = get_pod_vecs(poddim)
        lyitVy = facmy.solve_Ft(ypodvecs)
        red_realize_sol, red_realize_output, red_probfems, red_plotit \
            = get_red_problem(lyitVy)
        red_cmat = red_probfems['cmat']

        if checkredmod:
            nulist = [basenu]*uncdims
            redv = red_realize_sol(nulist)
            red_plotit(vvec=redv, pvdfile=redsolfile, plotplease=plotplease)
            print('N{1}pod{2}red_y(basenu)={0}'.format(red_cmat.dot(redv),
                                                       meshlevel, poddim))

        if pcepod:
            dimsrpcexpyl = []
            for pcedim in pcedimlist:
                abscissae, weights, compredexpv = mpu.\
                    setup_pce(distribution='uniform',
                              distrpars=dict(a=nua, b=nub),
                              pcedim=pcedim, uncdims=uncdims)
                redysoltens = mpu.\
                    run_pce_sim_separable(solfunc=red_realize_output,
                                          uncdims=uncdims, abscissae=abscissae)
                redpcexpy = compredexpv(redysoltens)
                dimsrpcexpyl.append(redpcexpy)
                print('pcedim={0:2.0f}, poddim={2:2.0f}, exypce={1}'.
                      format(pcedim, redpcexpy-pcexpy, poddim))
            rpcesxpyl.append(dimsrpcexpyl)
        if mcpod:
            varinu = basenu+(varib-varia)*np.random.rand(redmcruns, uncdims)
            expvnu = np.average(varinu, axis=0)
            print('expected value of nu: ', expvnu)
            varinulist = varinu.tolist()
            mcout, rmcxpy, expvnu = mpu.\
                run_mc_sim(varinulist, red_realize_output)
            cmpwmc += 'poddim={0:2.0f}, rmcxpy-mcxpy={1}'.\
                format(poddim, rmcxpy-mcxpy)
            try:
                cmpwpce += 'poddim={0:2.0f}, rmcxpy-pcexpy={1}'.\
                    format(poddim, rmcxpy-pcexpy)
            except NameError:
                pass

            print('mcruns={0:2.0f}, poddim={2:2.0f}, rmcxpy-mcxpy={1}'.
                  format(redmcruns, rmcxpy-mcxpy, poddim))
            rmcxpyl.append(rmcxpy)

    plt.show()


if __name__ == '__main__':
    problem = 'cylinder'
    meshlevel = 7
    mcruns = 10  # 200
    pcedimlist = [3]  # , 3, 4, 5]  # , 7]
    mcplease = False
    pceplease = False
    plotplease = False
    mcpod = False
    pcepod = False
    # ## make it come true
    # mcplease = True
    # pceplease = True
    plotplease = True
    # pcepod = True
    # mcpod = True
    basisfrom = 'mc'
    basisfrom = 'pce'
    simit(mcruns=mcruns, pcedimlist=pcedimlist, problem=problem,
          meshlevel=meshlevel,
          plotplease=plotplease, basisfrom=basisfrom,
          mcplease=mcplease, pceplease=pceplease, mcpod=mcpod, pcepod=pcepod)
