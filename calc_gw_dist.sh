#!/bin/bash
#PJM -L rscgrp=medium-a
#PJM -L node=1
#PJM --mpi proc=576
#PJM -L elapse=0:01:00
#PJM -g gs58
#PJM -o /work/gs58/s58007/rna_GWclustering/job/log
#PJM -e /work/gs58/s58007/rna_GWclustering/job/err
#PJM -s 
. ~/.bashrc
module load odyssey

python="/work/gs58/s58007/app/anaconda3/envs/GW_on_PDB/bin/python"
script="/work/gs58/s58007/rna_GWclustering/compute_gw.py"

# Calculate the pair index from PJM_ARRAY_INDEX
PAIR_INDEX=${PJM_ARRAY_INDEX}
MATRIX_FILES=(/work/gs58/s58007/rna_GWclustering/data/internal_distmat/*csv)
TOTAL_FILES=${#MATRIX_FILES[@]}

# Calculate the indices i and j for the pair (MATRIX_FILES[i], MATRIX_FILES[j])
I=$(awk -v idx="$PAIR_INDEX" 'BEGIN{
    for (i = 0; i < '"$TOTAL_FILES"'; i++) {
        for (j = i + 1; j < '"$TOTAL_FILES"'; j++) {
            if (count == idx) {
                print i;
                exit;
            }
            count++;
        }
    }
}')

J=$(awk -v idx="$PAIR_INDEX" -v i="$I" 'BEGIN{
    count=0;
    for (j = i + 1; j < '"$TOTAL_FILES"'; j++) {
        if (count == idx) {
            print j;
            exit;
        }
        count++;
    }
}')



# Execute the Python script
# $python /work/gs58/s58007/rna_GWclustering/compute_gw.py "${MATRIX_FILES[$I]}" "${MATRIX_FILES[$J]}"