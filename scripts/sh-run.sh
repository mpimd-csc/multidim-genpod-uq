source addpypath.sh

MESH=4
MC=0
PCE=1
PCEDIMS='2-3'  # -4-5'
PODDIMS='3-6-9-12-15-16'
PCEPOD=1
MCPOD=0
BASISFROM=pce
BASISFROM=mc
BASISFROM=rb
MCRUNS=1000
PCESNAPDIM=2
NPCESNAP=$(($PCESNAPDIM**4))
MCSNAP=$((5*1*$NPCESNAP))
RBSNAP=$(($PCESNAPDIM**4))
VARINU='[3,7]e-4'
DIST='uniform'
# DIST='beta-2-5'
# PCEXPY=0.88102114
# value of PCE(5) for MESH=10 and VARINU='[3,7]e-4'
NPROCS=4
TIMINGS=1
LOGFILE=alldump
echo 'tail -f logs/'$LOGFILE

python3 run_the_sims.py --mesh $MESH \
    --mc $MC --mcruns $MCRUNS \
    --pce $PCE --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --poddims $PODDIMS --rombase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM \
    --varinu $VARINU --distribution $DIST \
    --rbsnap $RBSNAP \
    --mcpod $MCPOD --mcsnap $MCSNAP
    # >> logs/$LOGFILE
    # --mcpod 1 \
    # --mc 1000 \
