# source setgmesh.sh
# tools for mesh generation
# GMSH=gmsh
# DOLFINCONVERT=dolfin-convert
alias gmsh='~/software/other/gmsh-git-Linux64/bin/gmsh'
alias dolfin-convert='~/software/gits/fenics-components/dolfin/scripts/dolfin-convert/dolfin-convert'

GZIP=gzip
PYTHON=python3

# refinement levels
SCALES=("2.0" "1.0")  # "0.75")  #  "0.5" "0.3" "0.2" "0.1")

# geo file
GEO_FILE=5-segs.geo

# prefix for mesh files
MSH_PREFIX=5-segs

## CLEAR FILES AND START MESH GENERATION
find . -name "*.xml" -exec rm  {} \;
find . -name "*.msh" -exec rm  {} \;
find . -name "*.xml.gz" -exec rm  {} \;
find . -name "*.pos" -exec rm  {} \;
find . -name "*.pvd" -exec rm  {} \;
find . -name "*.vtu" -exec rm  {} \;

# convert for level
LVL=1

# create meshes for different levels
for SCALE in "${SCALES[@]}"; do

    echo "DO MESHING FOR SEGMENTS LEVEL=$LVL and SCALE=$SCALE"

    # create mesh with gmsh, convert to xmkl and zip it
    gmsh ${GEO_FILE} -check -clscale ${SCALE} -2 -o ${MSH_PREFIX}_lvl$LVL.msh
    # ~/software/other/gmsh-git-Linux64/bin/gmsh karman2D-rotcyl-bm.geo -check -clscale .5 -2 -o karman2D-rotcyl-bm.msh

    # convert it to xml
	dolfin-convert ${MSH_PREFIX}_lvl$LVL.msh ${MSH_PREFIX}_lvl$LVL.xml

    # zip file
    ${GZIP} -f ${MSH_PREFIX}_lvl${LVL}.xml
    ${GZIP} -f ${MSH_PREFIX}_lvl${LVL}_facet_region.xml
    ${GZIP} -f ${MSH_PREFIX}_lvl${LVL}_physical_region.xml

    # convert to paraview
    ${PYTHON} mesh_to_paraview.py           ${MSH_PREFIX}_lvl${LVL}.xml.gz  ${MSH_PREFIX}_lvl${LVL}.pvd
    ${PYTHON} meshfunction_to_paraview.py   ${MSH_PREFIX}_lvl${LVL}.xml.gz  ${MSH_PREFIX}_lvl${LVL}_facet_region.xml.gz ${MSH_PREFIX}_lvl${LVL}_facet_region.pvd
    ${PYTHON} meshfunction_to_paraview.py   ${MSH_PREFIX}_lvl${LVL}.xml.gz  ${MSH_PREFIX}_lvl${LVL}_physical_region.xml.gz ${MSH_PREFIX}_lvl${LVL}_physical_region.pvd

    # increment level counter
    LVL=$[$LVL+1]

    echo "#####################################################################\n"

# an unrolled version
gmsh ${GEO_FILE} -0


done
