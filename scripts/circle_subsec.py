import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import matplotlib.pyplot as plt

import dolfin

import dolfin_navier_scipy.dolfin_to_sparrays as dts

Nrgs = 5
meshprfx = '../mesh/5-segs_lvl1'
meshfile = meshprfx + '.xml.gz'
physregs = meshprfx + '_physical_region.xml.gz'
fctsregs = meshprfx + '_facet_region.xml.gz'


def get_problem():
    mesh = dolfin.Mesh(meshfile)
    boundaries = dolfin.MeshFunction('size_t', mesh, fctsregs)
    subdomains = dolfin.MeshFunction('size_t', mesh, physregs)

    dx = dolfin.Measure('dx', subdomain_data=subdomains)

    V = dolfin.FunctionSpace(mesh, 'CG', 1)

    gzero = dolfin.Constant((0))
    zbcs = dolfin.DirichletBC(V, gzero, boundaries, Nrgs)
    bcinds, bcvals = dts.unroll_dlfn_dbcs([zbcs])

    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    # the mass matrix
    mmat = dolfin.assemble(dolfin.inner(u, v)*dolfin.dx)
    mmat = dts.mat_dolfin2sparse(mmat)
    mmat, _, bcinidx = dts.condense_velmatsbybcs(mmat, return_bcinfo=True,
                                                 dbcinds=bcinds,
                                                 dbcvals=bcvals)

    ininds = bcinidx['ininds']

    lplclst = []
    kk = 1
    nv = len(ininds)

    for kk in range(Nrgs):
        # assemble all mit nu=1
        akform = dolfin.assemble((1.*dolfin.inner(dolfin.grad(u),
                                 dolfin.grad(v)))*dx(kk))
        Akmat = dts.mat_dolfin2sparse(akform)
        Akmat.eliminate_zeros()
        Akmat, _ = dts.condense_velmatsbybcs(Akmat, invinds=ininds,
                                             dbcinds=bcinds, dbcvals=bcvals)
        lplclst.append(Akmat)

    def realize_linop(nulist):
        amat = sps.csr_matrix((nv, nv))
        for kk, knu in enumerate(nulist):
            amat = amat + knu*lplclst[kk]
        return amat

    examplerhsexp = dolfin.\
        Expression("sin(2*pi*(pow(x[0],2)+pow(x[1],2)))*sin(pi*4*x[0])",
                   degree=1)
    rhs = dolfin.assemble(v*examplerhsexp*dx(2))
    examplerhsvec = rhs.get_local()[ininds]

    def realize_sol(nulist, rhs=None):
        rhsvec = examplerhsvec if rhs is None else rhs
        amat = realize_linop(nulist)
        solvec = spsla.spsolve(amat, rhsvec).reshape((nv, 1))
        return solvec

    problemfems = dict(mmat=mmat,
                       bcinds=bcinds, bcvals=bcvals, ininds=ininds)

    def plotit(vvec=None, vfun=None, fignum=1):
        if vfun is None:
            vfun = dts.expand_dolfunc(vvec, bcinds=bcinds, bcvals=bcvals,
                                      ininds=ininds, V=V)
        plt.figure(fignum)
        dolfin.plot(vfun)
        plt.show()

    def realize_output(nulist, rhs=None, plotfignum=None):
        solvec = realize_sol(nulist, rhs=rhs)
        solv = dts.expand_dolfunc(solvec, bcinds=bcinds, bcvals=bcvals,
                                  ininds=ininds, V=V)
        if plotfignum is not None:
            plotit(vfun=solv, fignum=plotfignum)
        output = np.sqrt(dolfin.assemble(solv*solv*dx(0)))
        return output

    return realize_linop, realize_sol, realize_output, problemfems
