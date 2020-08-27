export PYTHONPATH="$PYTHONPATH:../"

# ## Set the parameters
VARINU='[3,7]e-4'
MESH=10

NPROCS=4  # number of parallel processes

# ## PCE simulation
PCE=1
PCEDIMS='2-3-4-5'
TIMINGS=5  # set this to 1 if only the results are wanted
python3 run_the_sims.py --mesh $MESH \
    --pce $PCE --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --varinu $VARINU \
    --plotpcepoddifddim $PCEPLOTDIFFDIM

# ## POD Approximation with the PCE basis
PCEDIMS='2-3-4-5'
PODDIMS='3-6-9-12-15-16'
PCEPOD=1
BASISFROM=pce
PCESNAPDIM=2
TIMINGS=5  # set this to 1 if only the results are wanted
PCEXPY=0.88102114  
# value of PCE(5) for MESH=10 and VARINU='[3,7]e-4'
# not needed but used to print the error directly to the screen
python3 run_the_sims.py --mesh $MESH \
    --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --poddims $PODDIMS --podbase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM --pcexpy $PCEXPY \
    --varinu $VARINU
# Postprocessing
python3 post-process.py N10nu3.00e-04--7.00e-04_pcepod_bfpce2.json

# ## POD Approximation with the random basis
PCEDIMS='2-3-4-5'
PODDIMS='3-6-9-12-15-16'
PCEPOD=1
BASISFROM=mc
PCESNAPDIM=2
NPCESNAP=$(($PCESNAPDIM**4))
MCSNAP=$((1*$NPCESNAP))
TIMINGS=5  # here the 5 is needed to compute the median
PCEXPY=0.88102114  
# value of PCE(5) for MESH=10 and VARINU='[3,7]e-4'
# not needed but used to print the error directly to the screen
python3 run_the_sims.py --mesh $MESH \
    --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --poddims $PODDIMS --podbase $BASISFROM \
    --mcsnap $MCSNAP \
    --pcepod $PCEPOD --pcexpy $PCEXPY \
    --varinu $VARINU
python3 post-process.py N10nu3.00e-04--7.00e-04_pcepod_bfmc16_runs5.json
# ## the results with 10 runs
TIMINGS=10
python3 run_the_sims.py --mesh $MESH \
    --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --poddims $PODDIMS --podbase $BASISFROM \
    --mcsnap $MCSNAP \
    --pcepod $PCEPOD --pcexpy $PCEXPY \
    --varinu $VARINU
python3 post-process.py N10nu3.00e-04--7.00e-04_pcepod_bfmc16_runs10.json
# ## the results with 80 snapshots
TIMINGS=5
MCSNAP=$((5*$NPCESNAP))
python3 run_the_sims.py --mesh $MESH \
    --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --poddims $PODDIMS --podbase $BASISFROM \
    --mcsnap $MCSNAP \
    --pcepod $PCEPOD --pcexpy $PCEXPY \
    --varinu $VARINU
python3 post-process.py N10nu3.00e-04--7.00e-04_pcepod_bfmc80_runs5.json
