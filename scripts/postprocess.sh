#!/bin/bash

# checking the PCE diffs and generating the plots
python3 pce-diffs.py --distribution 'uniform'
python3 pce-diffs.py --distribution 'beta-2-5'
python3 post-process.py --distribution 'uniform'
python3 post-process.py --distribution 'beta-2-5'


# the mesh test
python3 test_the_meshs.py --smlmesh 8 --lrgmesh 14
