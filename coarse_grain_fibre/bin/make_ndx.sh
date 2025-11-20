# Create initial index file
echo "q" | gmx make_ndx -f $1 -o tmp.ndx

# Create a temporary script to generate the selection commands
cat > make_selections.sh << 'EOF'
#!/bin/bash

# System details
atoms_per_mol=28
n_mol_per_block=50
total_molecules=200

# Create selections for iso-R and iso-S
start=$((0 * n_mol_per_block * atoms_per_mol + 1))
end=$(((0+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end | "
start=$((2 * n_mol_per_block * atoms_per_mol + 1))
end=$(((2+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end\n"

printf "name 5 iso-R\n"

start=$((1 * n_mol_per_block * atoms_per_mol + 1))
end=$(((1+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end | "
start=$((3 * n_mol_per_block * atoms_per_mol + 1))
end=$(((3+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end\n"

printf "name 6 iso-S\n"

# VMD to GROMACS index conversion (VMD is 1-indexed, GROMACS is 0-indexed). 
# So we subtract 1 from each VMD index
lleg_gmx=(14 18 19 20)
rleg_gmx=(15 17 21 22)
hed_gmx=(24 25 26 27)

# Create selections for LLEG, RLEG, HED, and CORE
printf "a "
for mol in $(seq 0 $((total_molecules-1))); do
  for atom in "${lleg_gmx[@]}"; do
    printf "$((mol * atoms_per_mol + atom + 1)) "
  done
done
printf "\n name 7 LLEG\n"

printf "a "
for mol in $(seq 0 $((total_molecules-1))); do
  for atom in "${rleg_gmx[@]}"; do
    printf "$((mol * atoms_per_mol + atom + 1)) "
  done
done
printf "\n name 8 RLEG\n"

printf "a "
for mol in $(seq 0 $((total_molecules-1))); do
  for atom in "${hed_gmx[@]}"; do
    printf "$((mol * atoms_per_mol + atom + 1)) "
  done
done
printf "\n name 9 HED\n"

# Generate CORE selection by excluding other groups
printf "a 1-$((total_molecules * atoms_per_mol)) & ! 7 & ! 8 & ! 9\n"
printf "name 10 CORE\n"
printf "del 2-4\n"

echo "q"
EOF

chmod +x make_selections.sh
# Generate the selection commands and pipe them to make_ndx
./make_selections.sh | gmx make_ndx -f $1 -n tmp.ndx -o $2