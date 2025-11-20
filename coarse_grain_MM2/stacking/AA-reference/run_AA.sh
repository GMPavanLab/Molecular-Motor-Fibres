#!/bin/bash


STRUCTURE=$1
TOPOLOGY=$2
NCORES=$3
GPU_ID=$4
mkdir -p tmp/
mkdir -p logs/
mkdir -p output/
MOL=MM2-AA

cat << EOF > plumed_MM2-AA.dat
WHOLEMOLECULES ENTITY0=1-49, ENTITY1=50-98

com1: CENTER ATOMS=1-49
com2: CENTER ATOMS=50-98

dist: DISTANCE ATOMS=com1,com2

metad: METAD ARG=dist PACE=500 HEIGHT=0.2 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 
    GRID_MIN=0.1 GRID_MAX=4.0 GRID_BIN=400 TEMP=300
    ...

PRINT STRIDE=500 ARG=dist,metad.bias FILE=COLVAR
EOF

cp $TOPOLOGY ${MOL}.top


gmx insert-molecules -ci $STRUCTURE -nmol 2 -box 4.0 4.0 4.0 -o tmp/${MOL}_box.gro
TOPOL=${MOL}.top
sed -i '/^MM2/s/1/2/' $TOPOL


gmx solvate -cp tmp/${MOL}_box.gro -cs spc216.gro -o tmp/${MOL}_solv.gro -p $TOPOL

mkdir -p em
echo "Energy minimisation"
gmx grompp -f mdp/emin.mdp -c tmp/${MOL}_solv.gro -p $TOPOL -o em/em_${MOL}.tpr -maxwarn 1
gmx mdrun -deffnm em/em_${MOL}  -nt $NCORES > logs/em_${MOL}_run.out 

printf '10\n0\n' | gmx energy -f em/em_${MOL}.edr -o output/em_${MOL}_potential.xvg
echo "Done, potential in logs"

mkdir -p nvt
echo "NVT"
gmx grompp -f mdp/nvt.mdp -c em/em_${MOL}.gro -r em/em_${MOL}.gro -p $TOPOL -o nvt/nvt_${MOL}.tpr -maxwarn 1
gmx mdrun -deffnm nvt/nvt_${MOL} -nt $NCORES -gpu_id $GPU_ID

printf '15\n0\n' | gmx energy -f nvt/nvt_${MOL}.edr -o output/nvt_${MOL}_temperature.xvg
echo "Done, temperature in logs"

mkdir -p npt
echo "NPT"
gmx grompp -f mdp/npt.mdp -c nvt/nvt_${MOL}.gro -r nvt/nvt_${MOL}.gro -t nvt/nvt_${MOL}.cpt -p $TOPOL -o npt/npt_${MOL}.tpr -maxwarn 1
gmx mdrun -deffnm npt/npt_${MOL} -nt $NCORES -gpu_id $GPU_ID

printf '17\n0\n' | gmx energy -f npt/npt_${MOL}.edr -o output/npt_${MOL}_pressure.xvg 
printf '23\n0\n' | gmx energy -f npt/npt_${MOL}.edr -o output/npt_${MOL}_density.xvg 

echo "Done, pressure and density in logs"

mkdir -p run
echo "MD Run 100 ns"
gmx grompp -f mdp/run.mdp -c npt/npt_${MOL}.gro -t npt/npt_${MOL}.cpt -p $TOPOL -o run/md_${MOL}.tpr -maxwarn 1

gmx mdrun -v -deffnm run/md_${MOL} -plumed plumed_MM2-AA.dat -gpu_id $GPU_ID -nt $NCORES

rm \#*
rm */\#*
rm bck*
plumed sum_hills --hills HILLS --stride 10000 --mintozero
echo "Done"
rm */\#*
