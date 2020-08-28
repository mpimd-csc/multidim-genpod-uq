# source setgmesh.sh
# tools for mesh generation
# GMSH=gmsh
# DOLFINCONVERT=dolfin-convert
# alias gmsh='~/software/other/gmsh-git-Linux64/bin/gmsh'
# alias dolfin-convert='~/software/gits/fenics-components/dolfin/scripts/dolfin-convert/dolfin-convert'

GZIP=gzip
PYTHON=python3

# refinement levels
SCALES=("2.0" "1.5" "1.2" "1.0" "0.9" "0.8" "0.72" "0.65" "0.59" "0.54" "0.5" "0.47" "0.45" "0.44")  #  "0.51" "0.5")
LVL=1

# geo file
GEO_FILE=4-segs-3D.geo

# prefix for mesh files
MSH_DIR=3D-mshs
MSH_PREFIX=4-segs-3D

## CLEAR FILES AND START MESH GENERATION
cd $MSH_DIR
# create meshes for different levels
for SCALE in "${SCALES[@]}"; do

    find . -name "*lvl${LVL}.xml" -exec rm  {} \;
    find . -name "*lvl${LVL}.msh" -exec rm  {} \;
    find . -name "*lvl${LVL}.xml.gz" -exec rm  {} \;
    find . -name "*lvl${LVL}.pos" -exec rm  {} \;
    find . -name "*lvl${LVL}.pvd" -exec rm  {} \;
    find . -name "*lvl${LVL}.vtu" -exec rm  {} \;


    echo "DO MESHING FOR SEGMENTS LEVEL=$LVL and SCALE=$SCALE"

    # create mesh with gmsh, convert to xmkl and zip it
    gmsh ../${GEO_FILE} -check -clscale ${SCALE} -3 -o ${MSH_PREFIX}_lvl$LVL.msh
    # ~/software/other/gmsh-git-Linux64/bin/gmsh karman2D-rotcyl-bm.geo -check -clscale .5 -2 -o karman2D-rotcyl-bm.msh

    # convert it to xml
	dolfin-convert ${MSH_PREFIX}_lvl$LVL.msh ${MSH_PREFIX}_lvl$LVL.xml

    # zip file
    ${GZIP} -f ${MSH_PREFIX}_lvl${LVL}.xml
    ${GZIP} -f ${MSH_PREFIX}_lvl${LVL}_facet_region.xml
    ${GZIP} -f ${MSH_PREFIX}_lvl${LVL}_physical_region.xml

    # convert to paraview
    ${PYTHON} ../mesh_to_paraview.py           ${MSH_PREFIX}_lvl${LVL}.xml.gz  ${MSH_PREFIX}_lvl${LVL}.pvd
    ${PYTHON} ../meshfunction_to_paraview.py   ${MSH_PREFIX}_lvl${LVL}.xml.gz  ${MSH_PREFIX}_lvl${LVL}_facet_region.xml.gz ${MSH_PREFIX}_lvl${LVL}_facet_region.pvd
    ${PYTHON} ../meshfunction_to_paraview.py   ${MSH_PREFIX}_lvl${LVL}.xml.gz  ${MSH_PREFIX}_lvl${LVL}_physical_region.xml.gz ${MSH_PREFIX}_lvl${LVL}_physical_region.pvd

    # ${PYTHON} trying_meshio.py
    # increment level counter
    LVL=$[$LVL+1]

    echo "#####################################################################\n"

    # an unrolled version
    gmsh ${GEO_FILE} -0

done
cd ..
