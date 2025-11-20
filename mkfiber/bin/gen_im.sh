#!/bin/bash
#UTILITY SCRIPT to generate xvg-files
mkdir -p ../output/NPT

printf '17\n0\n' | gmx energy -f $1 -o ../output/NPT/${2}_pressure.xvg
printf '22\n0\n' | gmx energy -f $1 -o ../output/NPT/${2}_volume.xvg
printf '21\n0\n' | gmx energy -f $1 -o ../output/NPT/${2}_z.xvg
printf '23\n0\n' | gmx energy -f $1 -o ../output/NPT/${2}_density.xvg
