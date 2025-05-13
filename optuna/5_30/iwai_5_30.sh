#!/bin/bash
#PJM -L rscgrp=a-pj24001724
#PJM -L node=1
#PJM --mpi proc=120
#PJM -L elapse=1:00:00
#PJM -j

module load intel
module load impi
module load vasp

#mpiexec ~/vasp_6.4.3_vtst_genkai_0725/bin/vasp_std >& log
#vef.pl
echo "0   0.047   $(awk -v min=-1000 -v max=-500 'BEGIN{srand(); print min+rand()*(max-min)}')   0" > fe.dat

WORK_DIR=$(grep '^WORK_DIR=' .env | cut -d '=' -f2)
NUM_CANDIDATES=$(grep '^NUM_CANDIDATES=' .env | cut -d '=' -f2 | tr -d '[:space:]')
if ! [[ "$NUM_CANDIDATES" =~ ^[0-9]+$ ]]; then
    echo "Error: NUM_CANDIDATES is not a valid integer." >&2
    exit 1
fi

read -r col1 col2 col3 <<< "$(awk 'END {print $1, $2, $3}' fe.dat)"

dir_name=$(basename "$PWD")
if [[ $dir_name =~ ^([0-9]+)_([0-9]+)$ ]]; then
    trial=${BASH_REMATCH[1]}
    ID=${BASH_REMATCH[2]}
else
    echo "Error: Directory name format is incorrect." >&2
    exit 1
fi

if (( $(echo "$col2 <= 0.05" | bc -l) )); then
    output="$trial $ID $col3"
else
    output="$trial $ID 1"
fi

echo "$output" >> "$WORK_DIR/job.log"

log_lines=$(wc -l < "$WORK_DIR/job.log")
if [ "$log_lines" -eq "$NUM_CANDIDATES" ]; then
    while IFS=' ' read -r trial ID energy; do
        awk -F, -v trial="$trial" -v ID="$ID" -v energy="$energy" 'BEGIN {OFS=","} {
            if ($1 == trial && $2 == ID) {
                $8 = energy
            }
            print
        }' "$WORK_DIR/data.csv" > "$WORK_DIR/data_temp.csv" && mv "$WORK_DIR/data_temp.csv" "$WORK_DIR/data.csv"
    done < "$WORK_DIR/job.log"
    > "$WORK_DIR/job.log"
    RSA_KEY_PATH="/home/pj24001724/ku40000345/rsa_key/busseiken_private"
    ssh -i "$RSA_KEY_PATH" ku40000345@genkai "cd $WORK_DIR; qsub run_optuna_py.sh"
fi
