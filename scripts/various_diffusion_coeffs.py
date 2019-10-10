import numpy as np
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
amat = sps.csr_matrix((nv, nv))
for kk in range(Nrgs):
    knu = nulist[kk]
    akform = dolfin.assemble((knu*dolfin.inner(dolfin.grad(u),
                             dolfin.grad(v)))*dx(kk))
    Akmat = dts.mat_dolfin2sparse(akform)
    Akmat.eliminate_zeros()
    Akmat, _ = dts.condense_velmatsbybcs(Akmat, invinds=ininds,
                                         dbcinds=bcinds, dbcvals=bcvals)
    print(kk, Akmat.nnz, np.linalg.norm(Akmat.data))
    lplclst.append(Akmat)
    amat = amat + Akmat

onedolfun = dolfin.Constant('1')
testplotfile = dolfin.File('checkit.pvd')
for kk in range(Nrgs):
    onedolfvec = dolfin.assemble(v*onedolfun*dx(kk))
    monevec = ((onedolfvec.get_local())[ininds].reshape((nv, 1)))
    onevec = spsla.spsolve(mmat, monevec).reshape((nv, 1))
    aone = lplclst[kk]*onevec
    print(kk, np.linalg.norm(aone))
    aonefun = dts.expand_vecnbc_dolfunc(V=V, vec=aone, invinds=ininds,
                                        bcindsl=[bcinds], bcvalsl=[bcvals])
    # dolfin.plot(aonefun)
    # plt.show()
    testplotfile << aonefun, kk
