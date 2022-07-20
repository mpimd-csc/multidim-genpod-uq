# ## Set the parameters
MESH=12
VARINU='[5,10]e-4'
NPROCS=2
PCEDIMS='2-3-4-5'
PODDIMS='1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16'

# ## toggle for the different regimes  ## #
DIST='beta-2-5'
# DIST='uniform'

# ## PCE simulation (also as benchmark)
PCE=1

python3 run_the_sims.py --mesh $MESH \
    --varinu $VARINU --nprocs $NPROCS \
    --dist $DIST \
    --pce $PCE --pcedims $PCEDIMS


# ## POD check (based on the PCE-based basis)
PCEPOD=1
BASISFROM=pce
PCESNAPDIM=2

python3 run_the_sims.py --mesh $MESH \
    --varinu $VARINU \
    --distribution $DIST \
    --nprocs $NPROCS --timings $TIMINGS \
    --pcedims $PCEDIMS \
    --poddims $PODDIMS --rombase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM

# ## POD check (based on the random-snapshot basis)
MCSNAP=$((1*2**4))
PCEPOD=1
BASISFROM=mc

python3 run_the_sims.py --mesh $MESH \
    --varinu $VARINU \
    --distribution $DIST \
    --nprocs $NPROCS --timings $TIMINGS \
    --pcedims $PCEDIMS \
    --poddims $PODDIMS --rombase $BASISFROM \
    --mcsnap $MCSNAP \
    --pcepod $PCEPOD

# ## POD check (based on the reduced-basis basis)
RBSNAP=$((1*2**4))
PCEPOD=1
BASISFROM=rb

python3 run_the_sims.py --mesh $MESH \
    --varinu $VARINU \
    --distribution $DIST \
    --nprocs $NPROCS --timings $TIMINGS \
    --pcedims $PCEDIMS \
    --poddims $PODDIMS --rombase $BASISFROM \
    --rbsnap $RBSNAP \
    --pcepod $PCEPOD

# ## POD check (based on the random-snapshot basis -- more snapshots)
MCSNAP=$((1*2**4))
PCEPOD=1
BASISFROM=mc

python3 run_the_sims.py --mesh $MESH \
    --varinu $VARINU \
    --distribution $DIST \
    --nprocs $NPROCS --timings $TIMINGS \
    --pcedims $PCEDIMS \
    --poddims $PODDIMS --rombase $BASISFROM \
    --mcsnap $MCSNAP \
    --pcepod $PCEPOD

# ## POD check (based on the reduced-basis basis -- larger training set)
RBSNAP=$((1*2**4))
PCEPOD=1
BASISFROM=rb

python3 run_the_sims.py --mesh $MESH \
    --varinu $VARINU \
    --distribution $DIST \
    --nprocs $NPROCS --timings $TIMINGS \
    --pcedims $PCEDIMS \
    --poddims $PODDIMS --rombase $BASISFROM \
    --rbsnap $RBSNAP \
    --pcepod $PCEPOD
