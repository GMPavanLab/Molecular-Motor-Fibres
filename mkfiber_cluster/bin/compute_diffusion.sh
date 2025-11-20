#!/bin/bash

gro=$1
tpr=$2
xtc=$3
outname=$4
NATOMS=28
# Calculate atom range
echo $outname
# atoms per fibre = 220 * NATOMS
start1=$((220 * $5 * NATOMS + 1))
end1=$((220 * $6 * NATOMS))

mkdir -p ../tmp

if [ -n "$7" ] && [ -n "$8" ]; then
    # --- two ranges case: [start1,end1] U [start2,end2] ---
    start2=$((220 * $7 * NATOMS + 1))
    end2=$((220 * $8 * NATOMS))
    gmx make_ndx -f "$gro" -o ../tmp/${outname}.ndx << EOF
a ${start1}-${end1}
a ${start2}-${end2}
8 | 9
del 8
del 8
q
EOF

else
    # --- single range case: original behaviour ---
    gmx make_ndx -f "$gro" -o ../tmp/${outname}.ndx << EOF
a ${start1}-${end1}
q
EOF
fi
# Extract fiber with PBC mol
printf "8\n" | gmx trjconv -f $xtc -s $tpr -n ../tmp/${outname}.ndx -o ../tmp/${outname}_raw.xtc -pbc mol -dt 1000
echo 8 | gmx trjconv -f $gro -s $tpr -n ../tmp/${outname}.ndx -o ../tmp/${outname}_raw.gro -e 0 -pbc mol

# Apply nojump
echo 0 | gmx trjconv -f ../tmp/${outname}_raw.xtc -s ../tmp/${outname}_raw.gro -o ../tmp/${outname}_nojump.xtc -pbc nojump

# Run diffusion analysis
python compute_diffusion.py ../tmp/${outname}_raw.gro ../tmp/${outname}_nojump.xtc --start-fit 0.1 --end-fit 0.9 --molecule-indices 0 1 2 3 4 5 10 15 20 --validate-drift --output ../figures/diffusion/${outname}
