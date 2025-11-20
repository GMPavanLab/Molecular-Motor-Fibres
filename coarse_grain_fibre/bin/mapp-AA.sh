#!/bin/bash
cd ../mapp-traj
NAME=$1
echo 13 | gmx trjconv -f ../AA-reference/AA-run${NAME}.xtc -o AA-traj_whole${NAME}.xtc -s ../AA-reference/AA-run${NAME}.tpr \
   -pbc whole -dt 100

mkdir -p itp
python ../bin/hack_topol.py ../AA-reference/itp/motor2m-R.itp itp/motor2m-R_hacked.itp 
python ../bin/hack_topol.py ../AA-reference/itp/motor2m-S.itp itp/motor2m-S_hacked.itp 
cp ../AA-reference/AA-run_hacked${NAME}.top topol_hacked${NAME}.top
cp ../AA-reference/itp/sol_ions_gaff.itp itp/

rm mapping_expanded.ndx
python ../bin/expand_mapping.py ../structure_files/mapping.ndx mapping_expanded.ndx 152 100

gmx grompp -p topol_hacked${NAME}.top -f ../AA-reference/AA-run.mdp -c ../AA-reference/AA-npt${NAME}.gro -o AA-run_hacked${NAME}.tpr -maxwarn 2
rm mdout.mdp # clean-up

no_of_beads=$(grep "\[" mapping_expanded.ndx | wc -l)
no_of_beads_minus_1=$( python -c "print( $no_of_beads - 1)" )

seq 0 $no_of_beads_minus_1 | gmx traj -f AA-traj_whole${NAME}.xtc -s AA-run_hacked${NAME}.tpr -oxt AA-mapped${NAME}.xtc \
  -n mapping_expanded.ndx  -ng ${no_of_beads} -com


seq 0 $no_of_beads_minus_1 | gmx traj -f AA-traj_whole${NAME}.xtc -s AA-run_hacked${NAME}.tpr \
                                      -n mapping_expanded.ndx -oxt CG-fibre${NAME}.gro -ng ${no_of_beads} -com -e 0

# Generate the CG tpr
gmx grompp -f mdp/CG_run.mdp -c CG-fibre${NAME}.gro -p ../structure_files/CG_system.top -o CG-run${NAME}.tpr -maxwarn 1

# Generate index file
bash ../bin/make_ndx.sh CG-fibre${NAME}.gro CG-fibre${NAME}.ndx

cd ../bin

rm *#*
