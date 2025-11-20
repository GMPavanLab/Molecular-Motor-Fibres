#!/bin/bash/
:'
Generates a 5x3 fibre cluster with infinite fibres along y-axis and infinite stacking along x-axis. Building the system is not simple and often causes instabillities. 
    Solved and partly equillibrated systems are stored in ex. system/5_layers/equilibrated.gro. These should more often be used as a starting point for simulations
    rather than this script.
'
GPU_ID=$1
NT=$2
NAME=5_layers
cp ../structure_files/CG_system.top ../tmp/system.top

bash packmol_structure.sh ../structure_files/fibre_structure.pdb ../tmp/system.gro \
    ../structure_files/CA.pdb ../tmp/system.top 5 $GPU_ID

bash build_system.sh -f ../tmp/system.gro -p ../tmp/system.top -deffnm $NAME -gpu_id $GPU_ID -nt $NT

bash equilibriate.sh -f ../system/${NAME}/compressed.gro -p ../system/${NAME}_topol.top -deffnm $NAME -gpu_id $GPU_ID -nt $NT

mkdir -p ../run/run
gmx grompp -f mdp/run.mdp -p ../system/${NAME}_topol_sol.top -c ../system/${NAME}/equilibrated.gro \
     -n ../run/${NAME}.ndx -o ../run/run/${NAME}.tpr
gmx mdrun -v -deffnm ../run/run/${NAME} -gpu_id $GPU_ID -nt $NT -pin on -pinoffset 48

bash compute_diffusion.sh ../run/run/${NAME}.gro ../run/run/${NAME}.tpr ../run/run/${NAME}.xtc ${NAME}_bundle 0 1
