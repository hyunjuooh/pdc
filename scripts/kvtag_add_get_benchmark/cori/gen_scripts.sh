#!/bin/bash
N_THREAD=NO
MAX_NODE=512
MAX_ATTR=1024
MAX_ATTRLEN=1000

for (( i = 1; i <= $MAX_NODE; i*=2 )); do
    mkdir -p $i
    for (( j = 1; j <= $MAX_ATTR; j*=4 )); do
        for (( k = 100; k <= $MAX_ATTRLEN; k*=10 )); do
            JOBNAME=kvtag_bench_${i}_${j}_${k}
            TARGET=./$i/$JOBNAME.sbatch
            cp template.sh $TARGET
            sed -i "s/JOBNAME/${JOBNAME}/g"           $TARGET
            sed -i "s/NODENUM/${i}/g"           $TARGET
            sed -i "s/ATTRNUM/${j}/g"           $TARGET
            sed -i "s/ATTRLEN/${k}/g"           $TARGET
            if [[ "$i" -gt "16" ]]; then
                sed -i "s/REG//g"           $TARGET
            else
                sed -i "s/DBG//g"           $TARGET
            fi
        done
    done
done