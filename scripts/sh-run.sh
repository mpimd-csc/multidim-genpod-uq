source addpypath.sh
MESH=14
MC=0
PCE=0
PCEDIMS='3-4-5'
PCEPOD=0
MCPOD=0
BASISFROM=pce
MCRUNS=10000
PCESNAPDIM=2
NPCESNAP=$(($PCESNAPDIM**5))
MCSNAP=$((5*$NPCESNAP))
LOGFILE=alldump
echo 'tail -f logs/'$LOGFILE

python3 run_the_sims.py --mesh $MESH \
    --mc $MC --mcruns $MCRUNS \
    --pce $PCE --pcedims $PCEDIMS \
    --podbase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM \
    --mcpod $MCPOD --mcsnap $MCSNAP # \
    >> logs/$LOGFILE
    # --mcpod 1 \
    # --mc 1000 \
(48:2A:E3:5B:08:CB)
