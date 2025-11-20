#!/bin/bash
: '
Creates AA-references for coarse graining by running a 1000 ns simulation.
INPUT:
    $1: int
        Number of cores to use
    $2: int
        GPU ID to use
    $3: str
        Name of the molecule to simulate. ../structure_files/${MOL}.(gro|top) must both exist
OUTPUT:
    -run/* : Gromacs output files for the simulation
'
NCORES=$1
GPU_ID=$2
mkdir -p tmp/
mkdir -p logs/
mkdir -p output/
MOL=$3 #MM2 or MM2b-R

cp ../structure_files/${MOL}.top ${MOL}.top
echo "RUNNING ${MOL}"
echo "Editing conf"
gmx editconf -f ../structure_files/${MOL}.gro -o tmp/${MOL}_box.gro -c -d 1.0 -bt dodecahedron

TOPOL=${MOL}.top

gmx solvate -cp tmp/${MOL}_box.gro -cs spc216.gro -o tmp/${MOL}_solv.gro -p $TOPOL

echo "Adding ions"
gmx grompp -f mdp/ions.mdp -c tmp/${MOL}_solv.gro -p $TOPOL -o tmp/ions.tpr
printf '6\n0\n' | gmx genion -s tmp/ions.tpr -o tmp/${MOL}_ions.gro -p $TOPOL -pname NA -nname CL -neutral

mkdir -p em
echo "Energy minimisation"
gmx grompp -f mdp/emin.mdp -c tmp/${MOL}_ions.gro -p $TOPOL -o em/em_${MOL}.tpr -maxwarn 1
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
echo "MD Run 1 ns"
gmx grompp -f mdp/run.mdp -c npt/npt_${MOL}.gro -t npt/npt_${MOL}.cpt -p $TOPOL -o run/md_${MOL}_run.tpr -maxwarn 1

gmx mdrun -s run/md_${MOL}_run.tpr -o run/md_${MOL}.cpt -c run/md_${MOL}.gro -x run/md_${MOL}.xtc -e run/md_${MOL}.edr -nt $NCORES -gpu_id $GPU_ID
echo "Done"
rm */\#*
