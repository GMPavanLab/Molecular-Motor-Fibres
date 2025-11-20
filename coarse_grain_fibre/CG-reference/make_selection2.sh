#!/bin/bash

# System details
atoms_per_mol=28
n_mol_per_block=50
total_molecules=200

printf "del 2\ndel 2\ndel 2\ndel 2\ndel 4\n"

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