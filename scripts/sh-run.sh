source addpypath.sh

MESH=3
MC=0
PCE=0
PCEDIMS='2-3-4-5'
PODDIMS='5-10-15-20'
PCEPOD=1
MCPOD=0
BASISFROM=pce
MCRUNS=1000
PCESNAPDIM=2
NPCESNAP=$(($PCESNAPDIM**4))
MCSNAP=$((5*$NPCESNAP))
VARINU='[3,7]e-3'
PCEXPY=0.88102114
# value of PCE(5) for MESH=10 and VARINU='[3,7]e-4'
NPROCS=4
TIMINGS=2
LOGFILE=alldump
echo 'tail -f logs/'$LOGFILE

python3 run_the_sims.py --mesh $MESH \
    --mc $MC --mcruns $MCRUNS \
    --pce $PCE --pcedims $PCEDIMS \
    --nprocs $NPROCS --timings $TIMINGS \
    --poddims $PODDIMS --podbase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM --pcexpy $PCEXPY \
    --varinu $VARINU \
    --mcpod $MCPOD --mcsnap $MCSNAP # \
    >> logs/$LOGFILE
    # --mcpod 1 \
    # --mc 1000 \
