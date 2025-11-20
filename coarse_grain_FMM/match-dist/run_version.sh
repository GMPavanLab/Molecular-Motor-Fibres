#!/bin/bash

version=$1
NCORES=$2
GPU_ID=$3
iso=R
solute=../mapp-traj/CG_motor2m-${iso}.gro
solvent_box=../structure_files/box_CG_W_eq.gro
solvent_name=W
solvent_atoms=1

cp ${version}/CG_system_base.top ${version}/CG_system.top
gmx solvate -cp ${solute} -cs ${solvent_box} -p ${version}/CG_system.top -o ${version}/initial.gro -box 6 6 4  # NEW!

echo "Adding ions"
gmx grompp -f mdp/ions.mdp -c ${version}/initial.gro -p ${version}/CG_system.top -o ${version}/ions.tpr -maxwarn 1
printf "6\n" | gmx genion -s ${version}/ions.tpr -o ${version}/initial.gro -p ${version}/CG_system.top -pname NA -neutral

INDEX=${version}/index.ndx
printf "2 | 3 | 4 | 5\nname 15 MOTOR\n6 | 7\nname 16 SOL_ION\nq\n" | gmx make_ndx -f ${version}/initial.gro -o ${INDEX}

mkdir -p ${version}/em
gmx grompp -p ${version}/CG_system.top -c ${version}/initial.gro -f mdp/em.mdp  -o ${version}/em/em.tpr -po ${version}/em/em.mdp  -maxwarn 1
gmx mdrun -v -deffnm ${version}/em/em -gpu_id $GPU_ID -nt $NCORES 

mkdir -p ${version}/output
printf "8\n0\n" | gmx energy -f ${version}/em/em.edr -o ${version}/output/em.xvg

mkdir -p ${version}/eq1
gmx grompp -p ${version}/CG_system.top -c ${version}/em/em.gro -n $INDEX -f mdp/eq1.mdp  -o ${version}/eq1/eq.tpr -po ${version}/eq1/eq.mdp  -maxwarn 1
gmx mdrun -deffnm ${version}/eq1/eq -gpu_id $GPU_ID -nt $NCORES 

for i in {2..5}; do
    mkdir -p ${version}/eq$i
    gmx grompp -p ${version}/CG_system.top -c ${version}/eq$(( $i - 1 ))/eq.gro -f mdp/eq${i}.mdp \
     -n $INDEX -o ${version}/eq${i}/eq.tpr -po ${version}/eq${i}/eq.mdp  -maxwarn 1
    gmx mdrun -v -deffnm ${version}/eq${i}/eq -gpu_id $GPU_ID -nt $NCORES 
    printf "11\n0\n" | gmx energy -f ${version}/eq${i}/eq.edr -o ${version}/output/eq${i}-temp.xvg
    printf "12\n0\n" | gmx energy -f ${version}/eq${i}/eq.edr -o ${version}/output/eq${i}-pressure.xvg
    printf "18\n0\n" | gmx energy -f ${version}/eq${i}/eq.edr -o ${version}/output/eq${i}-density.xvg
done

mkdir -p ${version}/run

gmx grompp -p ${version}/CG_system.top -n $INDEX -c ${version}/eq5/eq.gro -f mdp/run.mdp \
 -o ${version}/run/run.tpr -po ${version}/run/run.mdp 
gmx mdrun -v -deffnm ${version}/run/run -gpu_id $GPU_ID -nt $NCORES 

rm ${version}/*/\#*
rm ${version}/\#*
rm \#*