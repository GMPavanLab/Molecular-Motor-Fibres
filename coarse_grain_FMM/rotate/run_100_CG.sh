#!/bin/bash

echo "100 runs of CG" > converged.out
for i in {1..220}; do     
bash rotate_CG.sh 16 3 65 5000
A=$(tail -n 1 CG/COLVAR | awk '{print $8}')
I=$(tail -n 1 CG/COLVAR | awk '{print $9}')
B=$(tail -n 1 CG/COLVAR | awk '{print $10}')
CHECK=STUCK
if (( $(echo "($A > $B) && ($I > $B)" | bc -l) )); then
    CHECK=CONVERGED
fi
echo $CHECK >> converged.out
echo "######################################################################################"; 
echo $i $CHECK; 
echo "######################################################################################"; 

done