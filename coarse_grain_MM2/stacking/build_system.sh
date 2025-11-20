#!/bin/bash

solute=$1
solute_top=$2
solvent=$3

#Add two molecules and solve
echo 0 | gmx editconf -f $solute -o solute.gro -princ
cp $solute_top topol.top

gmx insert-molecules -ci $solute -nmol 2 -box 4 4 4 -o box.gro
sed -i '/^MM2/s/1/2/' topol.top
gmx solvate -cp box.gro -cs $solvent -p topol.top -o MM2_solv.gro
rm \#*
