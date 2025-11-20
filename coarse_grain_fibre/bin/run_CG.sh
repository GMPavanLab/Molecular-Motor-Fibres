#!/bin/bash

mkdir -p ../CG-reference
cd ../CG-reference

NCORES=16
NAME=$1
GPU_ID=$2
box_width=12
NMOLS=200
NATOMS=28
solute=../mapp-traj/CG-fibre${NAME}.gro
solvent_box=../structure_files/box_CG_W_eq.gro
solvent_name=W
solvent_atoms=1

cp ../structure_files/CG_system_base.top CG_system${NAME}.top
c=$(tail -n 1 ../mapp-traj/CG-fibre${NAME}.gro | awk '{print $NF + 0.2}')
gmx editconf -f ../mapp-traj/CG-fibre${NAME}.gro -c -o box${NAME}.gro  -box $box_width $box_width $c #-angles 90 90 120 -bt triclinic
gmx solvate -cp box${NAME}.gro -cs ${solvent_box} -p CG_system${NAME}.top -o box${NAME}.gro  
python ../bin/remove_water_inside_fiber.py box${NAME}.gro box_solv${NAME}.gro CG_system${NAME}.top 2.5

echo "Adding ions"
gmx grompp -f mdp/ions.mdp -c box_solv${NAME}.gro -p CG_system${NAME}.top -o ions${NAME}.tpr -maxwarn 1
printf "6\n" | gmx genion -s ions${NAME}.tpr -o box_solv${NAME}.gro -p CG_system${NAME}.top -pname NA -neutral

INDEX=index${NAME}.ndx
echo "q" | gmx make_ndx -f box_solv${NAME}.gro -o tmp${NAME}.ndx

# Create a temporary script to generate the selection commands
cat > make_selections.sh << 'EOF'
#!/bin/bash

# System details
atoms_per_mol=28
n_mol_per_block=50
total_molecules=200

for i in $(seq 2 12); do
    printf "del 2\n"
done

# Create selections for iso-R and iso-S
start=$((0 * n_mol_per_block * atoms_per_mol + 1))
end=$(((0+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end | "
start=$((2 * n_mol_per_block * atoms_per_mol + 1))
end=$(((2+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end\n"

printf "name 4 iso-R\n"

start=$((1 * n_mol_per_block * atoms_per_mol + 1))
end=$(((1+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end | "
start=$((3 * n_mol_per_block * atoms_per_mol + 1))
end=$(((3+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end\n"

printf "name 5 iso-S\n"

# VMD to GROMACS index conversion (VMD is 1-indexed, GROMACS is 0-indexed)
# So we subtract 1 from each VMD index
lleg_gmx=(14 16 18 19 20)
rleg_gmx=(15 17 21 22 23)
hed_gmx=(24 25 26 27)

# Create selections for LLEG, RLEG, HED, and CORE
for atom in "${lleg_gmx[@]}"; do
    printf "a "
    for mol in $(seq 0 $((total_molecules-1))); do
        printf "$((mol * atoms_per_mol + atom + 1)) "
  done
  printf "\n"
done
printf "6 | 7 | 8 | 9 | 10\nname 11 LLEG\n"
printf "del 10\ndel 9\ndel 8\ndel 7\n del 6\n"

for atom in "${rleg_gmx[@]}"; do
    printf "a "
    for mol in $(seq 0 $((total_molecules-1))); do
        printf "$((mol * atoms_per_mol + atom + 1)) "
  done
  printf "\n"
done
printf "7 | 8 | 9 | 10 | 11\nname 12 RLEG\n"
printf "del 11\ndel 10\ndel 9\ndel 8\ndel 7\n"


printf "a "
for mol in $(seq 0 $((total_molecules-1))); do
  for atom in "${hed_gmx[@]}"; do
    printf "$((mol * atoms_per_mol + atom + 1)) "
  done
done
printf "\n name 8 HED\n"

# Generate CORE selection by excluding other groups
printf "0 & ! 2 & ! 3 & ! 7 & ! 8 & ! 6\n"
echo "name 9 CORE"

printf "2 | 3\nname 10 SOL_ION\n"
printf "4 | 5\nname 11 MOTOR\n"
echo "q"
EOF

chmod +x make_selections.sh
# Generate the selection commands and pipe them to make_ndx
./make_selections.sh | gmx make_ndx -f box_solv${NAME}.gro -n tmp${NAME}.ndx -o $INDEX
rm tmp${NAME}.ndx


mkdir -p em
gmx grompp -p CG_system${NAME}.top -c box_solv${NAME}.gro -f mdp/em.mdp  -o em/em${NAME}.tpr -po em/em.mdp  -maxwarn 1
gmx mdrun -v -deffnm em/em${NAME} -nt 12 -pin on #-gpu_id $GPU_ID -nt $NCORES 

mkdir -p output
printf "8\n0\n" | gmx energy -f em/em${NAME}.edr -o output/em${NAME}.xvg

mkdir -p eq1
gmx grompp -p CG_system${NAME}.top -c em/em${NAME}.gro -n $INDEX -f mdp/eq1.mdp  -o eq1/eq${NAME}.tpr -po eq1/eq.mdp  -maxwarn 1
gmx mdrun -deffnm eq1/eq${NAME} -nt 12 -pin on -pinoffset 16 #-gpu_id $GPU_ID -nt $NCORES

for i in {2..5}; do
    mkdir -p eq$i
    gmx grompp -p CG_system${NAME}.top -c eq$(( $i - 1 ))/eq${NAME}.gro -f mdp/eq${i}.mdp \
     -n $INDEX -o eq${i}/eq${NAME}.tpr -po eq${i}/eq.mdp  -maxwarn 1
    gmx mdrun -v -deffnm eq${i}/eq${NAME} -nt 12 -nb cpu -pinoffset 16 # -gpu_id $GPU_ID -nt $NCORES 

done

mkdir -p run

gmx grompp -p CG_system${NAME}.top -n $INDEX -c eq5/eq${NAME}.gro -f mdp/run.mdp \
 -o run/run${NAME}.tpr -po run/run.mdp 
gmx mdrun -v -deffnm run/run${NAME} -gpu_id $GPU_ID -nt $NCORES

rm */\#*
rm \#*
