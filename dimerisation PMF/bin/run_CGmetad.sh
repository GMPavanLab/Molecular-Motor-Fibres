#!/bin/bash
: '
Runs a 500 ns metadynamics simulation on the COM between two CG molecules
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
N=14 # Number of atoms in each motor

mkdir -p ../metad
mkdir -p ../metad/COLVARS_${NAME}
mkdir -p ../tmp/
mkdir -p ../logs/
cp $TOPOL ../metad/${NAME}.top
gmx insert-molecules -f $STRUCTURE1 -ci $STRUCTURE2 -nmol 1 -box 3.5 3.5 4.0 -o ../tmp/${NAME}_box.gro
TOPOL=../metad/${NAME}.top
printf "${SCND_MOL}     1\n" >> $TOPOL

gmx solvate -cp ../tmp/${NAME}_box.gro -cs ../structure_files/box_CG_W_eq.gro -o ../tmp/${NAME}_solv.gro -p $TOPOL

mkdir -p ../metad/em
echo "Energy minimisation"
gmx grompp -p  $TOPOL -c ../tmp/${NAME}_solv.gro -f mdp-CG/em.mdp  -o ../metad/em/em_${NAME}.tpr -maxwarn 1
gmx mdrun -deffnm ../metad/em/em_${NAME} -gpu_id $GPU_ID -nt $NCORES

printf "8\n0\n" | gmx energy -f ../metad/em/em_${NAME}.edr -o ../logs/em_${NAME}_potential.xvg

mkdir -p ../metad/eq1
gmx grompp -p $TOPOL -c ../metad/em/em_${NAME}.gro -f mdp-CG/eq1.mdp  -o ../metad/eq1/eq_${NAME}.tpr  -maxwarn 1
gmx mdrun -deffnm ../metad/eq1/eq_${NAME} -gpu_id $GPU_ID -nt $NCORES

for i in {2..5}; do
    mkdir -p ../metad/eq$i
    gmx grompp -p $TOPOL -c ../metad/eq$(( $i - 1 ))/eq_${NAME}.gro -r ../metad/eq$(( $i - 1 ))/eq_${NAME}.gro -f mdp-CG/eq${i}.mdp \
     -o ../metad/eq${i}/eq_${NAME}.tpr  -maxwarn 1
    gmx mdrun -deffnm ../metad/eq${i}/eq_${NAME} -gpu_id $GPU_ID -nt $NCORES
    printf "12\n0\n" | gmx energy -f ../metad/eq${i}/eq_${NAME}.edr -o ../logs/eq${i}_${NAME}-temp.xvg
    printf "13\n0\n" | gmx energy -f ../metad/eq${i}/eq_${NAME}.edr -o ../logs/eq${i}_${NAME}-pressure.xvg
    printf "19\n0\n" | gmx energy -f ../metad/eq${i}/eq_${NAME}.edr -o ../logs/eq${i}_${NAME}-density.xvg
done

mkdir -p ../metad/run

cat << EOF > plumed_MM2-CG.dat
WHOLEMOLECULES ENTITY0=1-${N}, ENTITY1=$(($N + 1))-$((2 * $N)))

com1: CENTER ATOMS=1-${N}
com2: CENTER ATOMS=$(($N + 1))-$((2 * $N)))

dist: DISTANCE ATOMS=com1,com2

metad: METAD ARG=dist PACE=500 HEIGHT=0.2 SIGMA=0.05 ...
    FILE=../metad/COLVARS_${NAME}/HILLS BIASFACTOR=10 
    GRID_MIN=0.1 GRID_MAX=4.0 GRID_BIN=400 TEMP=300
    ...

PRINT STRIDE=500 ARG=dist,metad.bias FILE=../metad/COLVARS_${NAME}/COLVAR_raw_${NAME}
EOF

gmx grompp -p $TOPOL -c ../metad/eq5/eq_$NAME.gro -f mdp-CG/run.mdp \
 -o ../metad/run/md_${NAME}.tpr
gmx mdrun -v -deffnm ../metad/run/md_${NAME} -plumed plumed_MM2-CG.dat -gpu_id $GPU_ID -nt $NCORES

rm \#*
rm ../metad/*/\#*
rm ../metad/*/bck*
