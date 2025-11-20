#!/bin/bash

molecule=MOTOR
version=$1
size_1=${#molecule}
res=MOTOR

if [[ $version == "AA" ]]; then
    echo "1\n" | gmx gyrate -f ../AA-reference/run/md2m-R.xtc -s ../AA-reference/run/md2m-R_run.tpr \
    -n ../structure_files/index-R.ndx -b 100000 -e 200000 \
    -o gyrate-AA.xvg
else
    echo "15\n" | gmx gyrate -f ../match-dist/${version}/run/run.xtc -s ../match-dist/${version}/run/run.tpr \
    -n ../match-dist/${version}/index.ndx \
    -o gyrate-CG.xvg
fi