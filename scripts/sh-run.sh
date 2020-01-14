source addpypath.sh

MESH=10
MC=1
PCE=0
PCEDIMS='2-3-4-5'
PODDIMS='5-10-15-20'
PCEPOD=0
MCPOD=0
BASISFROM=pce
MCRUNS=100
PCESNAPDIM=2
NPCESNAP=$(($PCESNAPDIM**5))
MCSNAP=$((5*$NPCESNAP))
VARINU='[3,7]e-4'
NPROCS=6
LOGFILE=alldump
echo 'tail -f logs/'$LOGFILE

python3 run_the_sims.py --mesh $MESH \
    --mc $MC --mcruns $MCRUNS \
    --pce $PCE --pcedims $PCEDIMS \
    --nprocs $NPROCS \
    --poddims $PODDIMS --podbase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM \
    --varinu $VARINU \
    --mcpod $MCPOD --mcsnap $MCSNAP # \
    >> logs/$LOGFILE
    # --mcpod 1 \
    # --mc 1000 \
