source addpypath.sh
MESH=1
MC=0
PCE=0
PCEPOD=0
BASISFROM=mc
MCRUNS=10000
PCESNAPDIM=3
NPCESNAP=$(($PCESNAPDIM**5))
MCSNAP=$((5*$NPCESNAP))
echo $MCSNAP

python3 run_the_sims.py --mesh $MESH --mc $MC --pce $PCE --podbase $BASISFROM \
    --pcepod $PCEPOD --pcesnapdim $PCESNAPDIM --mcsnap $MCSNAP --mcruns $MCRUNS
    >> logs/alldump
    # --mcpod 1 \
    # --mc 1000 \
