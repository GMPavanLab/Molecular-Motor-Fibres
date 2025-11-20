#!/bin/bash
: '
Mapps the AA-reference to the CG model based on the mapping in ../structure_files/mapping_${MOL}.ndx. 
INPUT: 
    - $1 : MM2 | MM2b-R
      The molecule to be mapped
OUTPUTS:
    - mapped-${MOL}.xtc : Mapped trajectory of the AA-reference
    - CG-${MOL}.tpr : Gromacs input file for the mapped trajectory
    - CG-${MOL}.gro : Structure file for the CG model
'
MOL=$1 #MM2 or MM2b-R
TMP=MM2 # The version used in production used the same mappinng for the actuated motor as the ground state. Only changes the numbering of the atoms
traj=md_${MOL}.xtc
echo 0 | gmx trjconv -f ../AA-reference/run/${traj} -o AA-traj-${MOL}.xtc -s ../AA-reference/run/md_${MOL}_run.tpr -pbc whole

python hack_topol.py ../AA-reference/${MOL}.top ${MOL}_hacked.top

gmx grompp -p ${MOL}_hacked.top -f ../AA-reference/mdp/run.mdp -c ../AA-reference/npt/npt_${MOL}.gro -o AA-COG-${MOL}.tpr -maxwarn 2
rm mdout.mdp # clean-up

no_of_beads=$(grep "\[" ../structure_files/mapping_${TMP}.ndx | wc -l)
no_of_beads_minus_1=$( python -c "print( $no_of_beads - 1)" )

seq 0 $no_of_beads_minus_1 | gmx traj -f AA-traj-${MOL}.xtc -s AA-COG-${MOL}.tpr -oxt mapped-${MOL}.xtc \
  -n ../structure_files/mapping_${TMP}.ndx  -ng ${no_of_beads} -com


seq 0 $no_of_beads_minus_1 | gmx traj -f AA-traj-${MOL}.xtc -s AA-COG-${MOL}.tpr \
                                      -n ../structure_files/mapping_${TMP}.ndx -oxt CG-${MOL}.gro -ng ${no_of_beads} -com -e 0

# Generate the CG tpr
gmx grompp -f mdp/CG_run.mdp -c CG-${MOL}.gro -p ../structure_files/CG_system.top -o CG-${MOL}.tpr -maxwarn 1



rm \#*
