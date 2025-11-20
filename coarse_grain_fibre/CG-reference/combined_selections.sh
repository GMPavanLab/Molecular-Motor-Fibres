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
printf "name 9 iso-R\n"

start=$((1 * n_mol_per_block * atoms_per_mol + 1))
end=$(((1+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end | "
start=$((3 * n_mol_per_block * atoms_per_mol + 1))
end=$(((3+1) * n_mol_per_block * atoms_per_mol))
printf "a $start-$end\n"
printf "name 10 iso-S\n"

# Create SOL_ION group (assuming these are groups 0 and 1 after the deletions)
printf "5 | 6\n"
printf "name 11 SOL_ION\n"

printf "q\n"
