#!/bin/bash

for i in 1 2 3; do
    for scale in 9 10 11; do
        echo "scale ${scale}, i ${i}"
          bash make_fiber.sh -pdb ../../structure_files/pdb/motor2m-R.pdb ../../structure_files/pdb/motor2m-S.pdb\
            -top ../../structure_files/top/motor2m-R.top ../../structure_files/top/motor2m-S.top\
            --number_of_molecules 50 50 -deffnm scales_n${scale}_${i} -ns $scale -gpu_id 2 
    done
done