# import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import matplotlib.pyplot as plt

import dolfin

import dolfin_navier_scipy.dolfin_to_sparrays as dts

Nrgs = 5
nulist = [k for k in range(1, Nrgs+1)]

meshprfx = '../mesh/5-segs_lvl1'

meshfile = meshprfx + '.xml.gz'
physregs = meshprfx + '_physical_region.xml.gz'
fctsregs = meshprfx + '_facet_region.xml.gz'

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
                                             dbcinds=bcinds, dbcvals=bcvals)

# for kk in range(nrgs):
#     mmat = dolfin.assemble(dolfin.inner(u, v)*dx(kk))
#     mmat = dts.mat_dolfin2sparse(mmat)
#     mmat, _, bcinix = dts.condense_velmatsbybcs(mmat, return_bcinfo=True,
#                                                 dbcinds=bcinds,
#                                                 dbcvals=bcvals)
#     print(kk, mmat.nnz)
#     plt.figure(101)
#     plt.spy(M)
#     plt.show()
#     print(kk, M.nnz)

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


def _realize_amat(nulist):
    amat = sps.csr_matrix((nv, nv))
    for kk, knu in enumerate(nulist):
        amat = amat + knu*lplclst[kk]
    return amat

# onedolfun = dolfin.Constant('1')
# testplotfile = dolfin.File('checkit.pvd')
# for kk in range(Nrgs):
#     onedolfvec = dolfin.assemble(v*onedolfun*dx(kk))
#     monevec = ((onedolfvec.get_local())[ininds].reshape((nv, 1)))
#     onevec = spsla.spsolve(mmat, monevec).reshape((nv, 1))
#     aone = lplclst[kk]*onevec
#     print(kk, np.linalg.norm(aone))
#     aonefun = dts.expand_vecnbc_dolfunc(V=V, vec=aone, invinds=ininds,
#                                         bcindsl=[bcinds], bcvalsl=[bcvals])
#     # dolfin.plot(aonefun)
#     # plt.show()
#     testplotfile << aonefun, kk

# rhs = dolfin.Expression("sin(3.14159*(x[1]*x[1]+[x[1]*x[1]))", degree=1)
rhsexp = dolfin.\
    Expression("sin(2*pi*(pow(x[0],2)+pow(x[1],2)))*sin(pi*4*x[0])", degree=1)

rhs = dolfin.assemble(v*rhsexp*dolfin.dx)
rhsvec = rhs.get_local()[ininds]

amat = _realize_amat(nulist)
solvec = spsla.spsolve(amat, rhsvec).reshape((nv, 1))
solv = dts.expand_dolfunc(solvec, bcinds=bcinds, bcvals=bcvals,
                          ininds=ininds, V=V)

# rhs = dolfin.Expression("- 10*exp(- pow(x[1] - 0.5, 2))", degree=1)
# rhsfun = dolfin.interpolate(rhsexp, V)
dolfin.plot(solv)
plt.show()
