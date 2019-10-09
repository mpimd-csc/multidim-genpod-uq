r"""Convert a Fenics MeshFunction to Paraview File"""

from __future__ import print_function, division
import sys
from dolfin import Mesh, MeshFunction, File



if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: mesh_to_paraview.py MESH.xml MESHFUNCTION.xml MESHFUNCTION.pvd")
        sys.exit(2)

    meshfile = sys.argv[1]
    meshfunctionfile = sys.argv[2]
    pvdfile = sys.argv[3]

    mesh = Mesh(meshfile)
    meshfunction = MeshFunction("size_t", mesh, meshfunctionfile)
    pvd = File(pvdfile)

    pvd << meshfunction
