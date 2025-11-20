#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status
# Utility script to read trajectories from makefiber and analyse them

bash clean_trajectory.sh -s ../run/mdrun350_2m.pdb -f ../run/${1}.xtc -p ../run/topol.top -t ../run/${1}.tpr -o ${2}
python measure_parameters.py ../trajectories/${2}/${2}.pdb ../trajectories/${2}/${2}.xtc ${2} radius
python measure_parameters.py ../trajectories/${2}/${2}.pdb ../trajectories/${2}/${2}.xtc ${2} torsion