#!/bin/bash
: '
Runs a 500 ns metadynamics simulation on the COM between two AA molecules
INPUT:
    $1: Number of cores to run on 
    $2: GPU ID to run on
    $3: .gro - file
        Structure file for the first molecule
    $4: .gro -file
        Structure file for the second molecule
    $5: .top -file
        Topology file, with the first molecule in it
    $6: str
        Name of the second molecule, as named in the itp - files.
    $7: str
        Name of the simulation. Used for output control
OUTPUT:
    -metad/COLVAR_${NAME}/: COLVAR files for each simulation
    -metad/(em|eq*|npt|nvt|run): Gromacs output files for each simulation
'

NCORES=$1
GPU_ID=$2
STRUCTURE1=$3
STRUCTURE2=$4
TOPOL=$5
SCND_MOL=$6
NAME=$7

mkdir -p ../metad
mkdir -p ../metad/COLVARS_${NAME}
mkdir -p ../tmp/
mkdir -p ../logs/

cp ${TOPOL} ../metad/${NAME}.top

cat << EOF > plumed_MM2.dat
WHOLEMOLECULES ENTITY0=1-49, ENTITY1=50-98

com1: CENTER ATOMS=1-49
com2: CENTER ATOMS=50-98

dist: DISTANCE ATOMS=com1,com2

metad: METAD ARG=dist PACE=5000 HEIGHT=0.2 SIGMA=0.05 ...
    FILE=../metad/COLVARS_${NAME}/HILLS BIASFACTOR=10 
    GRID_MIN=0.1 GRID_MAX=4.0 GRID_BIN=400 TEMP=300
    ...

PRINT STRIDE=5000 ARG=dist,metad.bias FILE=../metad/COLVARS_${NAME}/COLVAR_raw_${NAME}
EOF

gmx insert-molecules -f $STRUCTURE1 -ci $STRUCTURE2 -nmol 1 -box 3.5 3.5 4.0 -o ../tmp/${NAME}_box.gro
TOPOL=../metad/${NAME}.top
printf "${SCND_MOL}     1\n" >> $TOPOL

gmx solvate -cp ../tmp/${NAME}_box.gro -cs spc216.gro -o ../tmp/${NAME}_solv.gro -p $TOPOL

mkdir -p ../metad/em
echo "Energy minimisation"
gmx grompp -f mdp-AA/emin.mdp -c ../tmp/${NAME}_solv.gro -p $TOPOL -o ../metad/em/em_${NAME}.tpr -maxwarn 1
gmx mdrun -deffnm ../metad/em/em_${NAME} -nt $NCORES -gpu_id $GPU_ID > ../logs/em_${NAME}_run.out 

printf '10\n0\n' | gmx energy -f ../metad/em/em_${NAME}.edr -o ../logs/em_${NAME}_potential.xvg
echo "Done, potential in logs"

mkdir -p ../metad/nvt
echo "NVT"
gmx grompp -f mdp-AA/nvt.mdp -c ../metad/em/em_${NAME}.gro -r ../metad/em/em_${NAME}.gro -p $TOPOL -o ../metad/nvt/nvt_${NAME}.tpr -maxwarn 1
gmx mdrun -deffnm ../metad/nvt/nvt_${NAME} -nt $NCORES -gpu_id $GPU_ID

printf '15\n0\n' | gmx energy -f ../metad/nvt/nvt_${NAME}.edr -o ../logs/nvt_${NAME}_temperature.xvg
echo "Done, temperature in logs"

mkdir -p ../metad/npt
echo "NPT"
gmx grompp -f mdp-AA/npt.mdp -c ../metad/nvt/nvt_${NAME}.gro -r ../metad/nvt/nvt_${NAME}.gro -t ../metad/nvt/nvt_${NAME}.cpt -p $TOPOL -o ../metad/npt/npt_${NAME}.tpr -maxwarn 1
gmx mdrun -deffnm ../metad/npt/npt_${NAME} -nt $NCORES -gpu_id $GPU_ID 

printf '17\n0\n' | gmx energy -f ../metad/npt/npt_${NAME}.edr -o ../logs/npt_${NAME}_pressure.xvg 
printf '23\n0\n' | gmx energy -f ../metad/npt/npt_${NAME}.edr -o ../logs/npt_${NAME}_density.xvg 

echo "Done, pressure and density in logs"

mkdir -p ../metad/run
echo "MD Run 500 ns"
gmx grompp -f mdp-AA/run.mdp -c ../metad/npt/npt_${NAME}.gro -t ../metad/npt/npt_${NAME}.cpt -p $TOPOL -o ../metad/run/md_${NAME}.tpr -maxwarn 1

gmx mdrun -deffnm ../metad/run/md_${NAME} -plumed plumed_MM2.dat -gpu_id $GPU_ID -nt $NCORES -pin on -pinoffset 32

rm ../tmp/\#*
rm ../metad/*/\#*
rm ../metad/bck*
echo "Done"
rm */\#*
