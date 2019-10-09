r"""Convert a Fenics Mesh to Paraview File"""

from __future__ import print_function, division
import sys
from dolfin import Mesh, File



if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: mesh_to_paraview.py MESH.xml MESH.pvd")
        sys.exit(2)

    meshfile = sys.argv[1]
    pvdfile = sys.argv[2]

    mesh = Mesh(meshfile)
    pvd = File(pvdfile)

    pvd << mesh
