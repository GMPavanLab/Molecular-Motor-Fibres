#!/bin/bash

system=$1
topol=$2
NCORES=$3
GPU_ID=$4

mkdir -p CG_metad
cp ${topol} CG_metad/CG_system.top
cd CG_metad

mkdir -p em
gmx grompp -p CG_system.top -c ../$system -f ../mdp/em.mdp  -o em/em.tpr -maxwarn 1
gmx mdrun -deffnm em/em -gpu_id $GPU_ID -nt $NCORES

mkdir -p output
printf "8\n0\n" | gmx energy -f em/em.edr -o output/em.xvg

mkdir -p eq1
gmx grompp -p CG_system.top -c em/em.gro -f ../mdp/eq1.mdp  -o eq1/eq.tpr  -maxwarn 1
gmx mdrun -deffnm eq1/eq -gpu_id $GPU_ID -nt $NCORES

for i in {2..5}; do
    mkdir -p eq$i
    gmx grompp -p CG_system.top -c eq$(( $i - 1 ))/eq.gro -r eq$(( $i - 1 ))/eq.gro -f ../mdp/eq${i}.mdp \
     -o eq${i}/eq.tpr -po eq${i}/eq.mdp  -maxwarn 1
    gmx mdrun -v -deffnm eq${i}/eq -gpu_id $GPU_ID -nt $NCORES
    printf "12\n0\n" | gmx energy -f eq${i}/eq.edr -o output/eq${i}-temp.xvg
    printf "13\n0\n" | gmx energy -f eq${i}/eq.edr -o output/eq${i}-pressure.xvg
    printf "19\n0\n" | gmx energy -f eq${i}/eq.edr -o output/eq${i}-density.xvg
done

mkdir -p run

cat << EOF > plumed_MM2-CG.dat
WHOLEMOLECULES ENTITY0=1-14, ENTITY1=15-28

com1: CENTER ATOMS=1-14
com2: CENTER ATOMS=15-28

dist: DISTANCE ATOMS=com1,com2

metad: METAD ARG=dist PACE=500 HEIGHT=0.2 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 
    GRID_MIN=0.1 GRID_MAX=4.0 GRID_BIN=400 TEMP=300
    ...

PRINT STRIDE=500 ARG=dist,metad.bias FILE=COLVAR
EOF

gmx grompp -p CG_system.top -c eq5/eq.gro -f ../mdp/run.mdp \
 -o run/run.tpr -po run/run.mdp 
gmx mdrun -v -deffnm run/run -plumed plumed_MM2-CG.dat -gpu_id $GPU_ID -nt $NCORES

rm \#*
rm */\#*
rm bck*
plumed sum_hills --hills HILLS --stride 10000 --mintozero
cd ..
