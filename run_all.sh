#!/bin/bash

MOLTYPE=2m #8 2m 11
GPU_ID=0

cd mkforce_field/bin
bash run_molecule.sh -s $MOLTYPE --head_charge 0 --leg_charge -1 --core_charge -1

cd ../../mkfiber/bin
bash run_all.sh -s $MOLTYPE -gpu_id $GPU_ID

cd ../../characterise_fibre/bin
bash run_all.sh
