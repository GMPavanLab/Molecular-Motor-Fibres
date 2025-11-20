#!/bin/bash

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <num_cores> <gpu_id> <tropic_label> [nrot]"
    exit 1
fi

NCORES=$1
GPU_ID=$2
TROPIC=$3
NROT=${4:-40}

bash run_AA.sh "$NCORES" "$GPU_ID" 150 "$TROPIC" "$NROT"

RUN_NAME="mdrun150_${TROPIC}"

bash clean_trajectory.sh \
    -s ../run/${RUN_NAME}.pdb \
    -f ../run/${RUN_NAME}.xtc \
    -p ../run/topol.top \
    -t ../run/${RUN_NAME}.tpr \
    -o "${TROPIC}"

python measure_parameters.py \
    ../trajectories/${TROPIC}/${TROPIC}.pdb \
    ../trajectories/${TROPIC}/${TROPIC}.xtc \
    "${TROPIC}" \
    radius
