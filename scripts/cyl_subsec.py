import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import numpy.linalg as npla
import matplotlib.pyplot as plt

import dolfin

import dolfin_navier_scipy.dolfin_to_sparrays as dts

dolfin.parameters['linear_algebra_backend'] = 'Eigen'

Nrgs = 4
meshprfx = '../mesh/3D-mshs/{0}-segs-3D_lvl1'.format(Nrgs)
meshfile = meshprfx + '.xml.gz'
physregs = meshprfx + '_physical_region.xml.gz'
fctsregs = meshprfx + '_facet_region.xml.gz'

# the physical entities of volumes, input faces, and output faces
volpes = [0 + k+1 for k in range(Nrgs)]
inppes = [10 + k+1 for k in range(Nrgs)]
outpes = [20 + k+1 for k in range(Nrgs)]


def get_problem():
    ''' the 3D cylinder
     * with Dirichlet control at the bottom
     * zero Neumann elsewhere
     * observation in a ring at the top
    '''

    mesh = dolfin.Mesh(meshfile)
    boundaries = dolfin.MeshFunction('size_t', mesh, fctsregs)
    subdomains = dolfin.MeshFunction('size_t', mesh, physregs)

    dx = dolfin.Measure('dx', subdomain_data=subdomains)
    ds = dolfin.Measure('ds', subdomain_data=boundaries)

    V = dolfin.FunctionSpace(mesh, 'CG', 1)

    # ## the boundary

    bcexp = dolfin.\
        Expression("sin(2*pi*(pow(x[0],2)+pow(x[1],2)))*sin(pi*4*x[0])",
                   degree=1)
    diribcs = []
    for pe in inppes:
        diribcs.append(dolfin.DirichletBC(V, bcexp, boundaries, pe))
    bcinds, bcvals = dts.unroll_dlfn_dbcs(diribcs)

    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    # the mass matrix
    mmat = dolfin.assemble(dolfin.inner(u, v)*dolfin.dx)
    mmat = dts.mat_dolfin2sparse(mmat)
    mmat, _, bcinidx = dts.condense_velmatsbybcs(mmat, return_bcinfo=True,
                                                 dbcinds=bcinds,
                                                 dbcvals=bcvals)
    ininds = bcinidx['ininds']

    lplclist, lplcrhslist = [], []
    kk = 1
    convfac = 1.

    convexp = dolfin.Expression(('(x[0]*x[0]+x[1]*x[1]-1)*x[1]',
                                 '(1-x[0]*x[0]-x[1]*x[1])*x[0]'), degree=1)
    convform = dolfin.assemble(v*dolfin.inner(convexp,
                               dolfin.grad(u))*dolfin.dx)
    convmat = convfac*dts.mat_dolfin2sparse(convform)
    convmat, convrhs = dts.\
        condense_velmatsbybcs(convmat, invinds=ininds,
                              dbcinds=bcinds, dbcvals=bcvals)

    for kk in volpes:
        # assemble all mit nu=1
        akform = dolfin.assemble((1.*dolfin.inner(dolfin.grad(u),
                                 dolfin.grad(v)))*dx(kk))
        Akmat = dts.mat_dolfin2sparse(akform)
        Akmat.eliminate_zeros()
        Akmat, krhs = dts.condense_velmatsbybcs(Akmat, invinds=ininds,
                                                dbcinds=bcinds, dbcvals=bcvals)
        lplclist.append(Akmat)
        lplcrhslist.appen(krhs)

    def realize_linop(nulist, lplclist=None, convmat=None):
        lnv = lplclist[0].shape[0]
        amat = sps.csr_matrix((lnv, lnv))
        for kk, knu in enumerate(nulist):
            amat = amat + knu*lplclist[kk]
        if convmat is None:
            return amat
        else:
            return amat + convmat

    def realize_rhs(nulist, lplcrhslist=None, convrhs=None):
        lnv = lplcrhslist[0].shape[0]
        rhs = np.zeros((lnv, 1))
        for kk, knu in enumerate(lplcrhslist):
            rhs = rhs + knu*lplclist[kk]
        if convmat is None:
            return rhs
        else:
            return rhs + convrhs

    def realize_sol(nulist, realize_amat=None, realize_rhs=None):
        amat = realize_amat(nulist)
        rhs = realize_rhs(nulist)
        if sps.issparse(amat):
            solvec = spsla.spsolve(amat, rhs).reshape((rhs.size, 1))
        else:
            solvec = npla.solve(amat, rhs).reshape((rhs.size, 1))
        return solvec

    def plotit(vvec=None, vfun=None, fignum=1):
        if vfun is None:
            vfun = dts.expand_dolfunc(vvec, bcinds=bcinds, bcvals=bcvals,
                                      ininds=ininds, V=V)
        plt.figure(fignum)
        dolfin.plot(vfun)
        plt.show()

    # ## Output

    obsoplist = []
    for pe in outpes:
        obsop = dolfin.assemble(u*ds(pe)).\
            get_local().reshape((1, V.dim()))
        obsopmat = sps.csc_matrix(obsop)[:, ininds]
        obsoplist.append(obsopmat/obsopmat.sum())

    cmat = sps.vstack(obsoplist)

    def realize_output(nulist, realize_sol=None, cmat=None,
                       plotfignum=None):
        solvec = realize_sol(nulist)
        if plotfignum is not None:
            solv = dts.expand_dolfunc(solvec, bcinds=bcinds, bcvals=bcvals,
                                      ininds=ininds, V=V)
            plotit(vfun=solv, fignum=plotfignum)
        output = cmat.dot(solvec)
        return output

    def full_realize_linop(nulist):
        return realize_linop(nulist, lplclist=lplclist, convmat=convmat)

    def full_realize_rhs(nulist):
        return realize_linop(nulist, lplcrhslist=lplcrhslist, convrhs=convrhs)

    def full_realize_sol(nulist):
        return realize_sol(nulist, realize_amat=full_realize_linop,
                           realize_rhs=full_realize_rhs)

    def full_realize_output(nulist, plotfignum=None):
        return realize_output(nulist, realize_sol=full_realize_sol, cmat=cmat,
                              plotfignum=plotfignum)

    problemfems = dict(mmat=mmat, cmat=cmat, realizeamat=full_realize_linop,
                       bcinds=bcinds, bcvals=bcvals, ininds=ininds,
                       realizerhs=full_realize_rhs)

    def get_red_prob(podmat):
        red_cmat = cmat.dot(podmat)

        red_lplclist = []
        red_lplcrhslist = []
        for kk in range(Nrgs):
            red_lplclist.append((podmat.T).dot(lplclist[kk].dot(podmat)))
            red_lplcrhslist.append((podmat.T).dot(lplcrhslist[kk]))

        red_convrhs = podmat.T.dot(convrhs)
        red_convmat = podmat.T.dot(convmat.dot(podmat))

        def red_realize_linop(nulist):
            return realize_linop(nulist, lplclist=red_lplclist,
                                 convmat=red_convmat)

        def red_realize_rhs(nulist):
            return realize_rhs(nulist, lplcrhslist=red_lplclist,
                               convrhs=red_convrhs)

        def red_realize_sol(nulist):
            return realize_sol(nulist, realize_amat=red_realize_linop,
                               rhs=red_realize_rhs)

        def red_realize_output(nulist, plotfignum=None):
            return realize_output(nulist, realize_sol=red_realize_sol,
                                  cmat=red_cmat, plotfignum=plotfignum)

        red_problemfems = dict(cmat=red_cmat, realizeamat=red_realize_linop,
                               realizerhs=red_realize_rhs,
                               convmat=red_convmat)
        return red_realize_sol, red_realize_output, red_problemfems

    return full_realize_sol, full_realize_output, problemfems, get_red_prob
