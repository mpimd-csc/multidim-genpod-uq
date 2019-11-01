import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import numpy.linalg as npla
import matplotlib.pyplot as plt

import dolfin

import dolfin_navier_scipy.dolfin_to_sparrays as dts

Nrgs = 5
meshprfx = '../mesh/5-segs_lvl1'
meshfile = meshprfx + '.xml.gz'
physregs = meshprfx + '_physical_region.xml.gz'
fctsregs = meshprfx + '_facet_region.xml.gz'


class CircObsvDomain(dolfin.SubDomain):

    def __init__(self, cx=0., cy=0., radius=.1):
        dolfin.SubDomain.__init__(self)
        self.cx = cx
        self.cy = cy
        self.rsqrd = radius**2

    def inside(self, x, on_boundary):
        xshftd = x[0]-self.cx
        yshftd = x[1]-self.cy
        insi = yshftd**2 + xshftd**2 < self.rsqrd
        # print(x, insi)
        return insi


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
    convfac = 1.

    convexp = dolfin.Expression(('(x[0]*x[0]+x[1]*x[1]-1)*x[1]',
                                 '(1-x[0]*x[0]-x[1]*x[1])*x[0]'), degree=1)
    convform = dolfin.assemble(v*dolfin.inner(convexp,
                               dolfin.grad(u))*dolfin.dx)
    convmat = convfac*dts.mat_dolfin2sparse(convform)
    convmat, _ = dts.condense_velmatsbybcs(convmat, invinds=ininds,
                                           dbcinds=bcinds, dbcvals=bcvals)

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
        return amat + convmat

    examplerhsexp = dolfin.\
        Expression("sin(2*pi*(pow(x[0],2)+pow(x[1],2)))*sin(pi*4*x[0])",
                   degree=1)
    rhs = dolfin.assemble(v*examplerhsexp*dx(0)+v*examplerhsexp*dx(2))
    examplerhsvec = rhs.get_local()[ininds]

    def realize_sol(nulist, rhs=None, podmat=None):
        if podmat is None:
            rhsvec = examplerhsvec if rhs is None else rhs
            amat = realize_linop(nulist)
            solvec = spsla.spsolve(amat, rhsvec).reshape((nv, 1))
            return solvec
        else:
            knv = podmat.shape[1]
            rhsvec = podmat.T.dot(examplerhsvec) if rhs is None else rhs
            amat = realize_linop(nulist)
            amat = podmat.T.dot(amat*podmat)
            solvec = npla.solve(amat, rhsvec).reshape((knv, 1))
            return podmat.dot(solvec)

    def plotit(vvec=None, vfun=None, fignum=1):
        if vfun is None:
            vfun = dts.expand_dolfunc(vvec, bcinds=bcinds, bcvals=bcvals,
                                      ininds=ininds, V=V)
        plt.figure(fignum)
        dolfin.plot(vfun)
        plt.show()

    # ## Output

    obsvr = .1
    obsarea = np.pi*obsvr**2
    obsvdom = CircObsvDomain(cx=.0, cy=0., radius=obsvr)
    cellmarkers = dolfin.MeshFunction('size_t', V.mesh(),
                                      V.mesh().topology().dim())
    obsvdom.mark(cellmarkers, 111)  # just 0 didnt work
    dxobs = dolfin.Measure('dx', subdomain_data=cellmarkers)
    cform = dolfin.assemble(1./obsarea*u*dxobs(111)).\
        get_local().reshape((1, V.dim()))
    cmat = sps.csc_matrix(cform)[:, ininds]

    def realize_output(nulist, rhs=None, plotfignum=None,
                       pointeva=False, podmat=None):
        solvec = realize_sol(nulist, rhs=rhs, podmat=podmat)
        if plotfignum is not None or pointeva:
            solv = dts.expand_dolfunc(solvec, bcinds=bcinds, bcvals=bcvals,
                                      ininds=ininds, V=V)
            if plotfignum is not None:
                plotit(vfun=solv, fignum=plotfignum)
                output = cmat.dot(solvec)
            if pointeva:
                midpt = dolfin.Point(0, 0)
                output = solv(midpt)
        else:
            output = cmat.dot(solvec)
            # output = np.sqrt(dolfin.assemble(solv*solv*dx(4)))
        return output

    problemfems = dict(mmat=mmat, cmat=cmat,
                       bcinds=bcinds, bcvals=bcvals, ininds=ininds)

    return realize_linop, realize_sol, realize_output, problemfems
