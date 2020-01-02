source addpypath.sh
MESH=1
MC=0
PCE=1
PCEPOD=1
BASISFROM=pce
MCRUNS=10000
PCESNAPDIM=2
NPCESNAP=$(($PCESNAPDIM**5))
MCSNAP=$((5*$NPCESNAP))
LOGFILE=alldump
echo 'tail -f logs/'$LOGFILE

python3 run_the_sims.py --mesh $MESH --mc $MC --pce $PCE --podbase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM --mcsnap $MCSNAP --mcruns $MCRUNS # \
    >> logs/$LOGFILE
    # --mcpod 1 \
    # --mc 1000 \
