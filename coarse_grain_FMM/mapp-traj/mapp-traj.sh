#!/bin/bash
: '
Mapps the AA-reference to the CG model based on the mapping in ../structure_files/mapping_${MOL}.ndx.
INPUT: 
    $1 : R | S
      The molecule isomer to be mapped
    Optional parameters:
    -xtc : trajectory file to be mapped
    -top : topology file to be mapped
    -tpr : input tpr file
    -mdp : input mdp file
    -start : input structure file for restarting the simulation
    -o|--output : output directory
OUTPUTS:
    - ${OUT}/mapped-${iso}.xtc : Mapped trajectory of the AA-reference
    - ${OUT}/CG-${iso}.tpr : Gromacs input file for the mapped trajectory
    - ${OUT}/CG_motor2m-${iso}.gro : Structure file for the CG model
'

iso=$1; shift
TOP=../AA-reference/topol2m-${iso}.top
TPR=../AA-reference/run/md2m-${iso}_run.tpr
TRAJ=../AA-reference/run/md2m-${iso}.xtc
MDP=../AA-reference/mdp/run.mdp
RESTART=../AA-reference/run/md2m-${iso}.gro
OUT=.
while [[ $# -gt 0 ]]; do
  case $1 in
    -xtc)
      shift #past argument
      TRAJ=$1
      shift #past value
      ;;
    -top)
      shift #past argument
      TOP=$1
      shift #past value
      ;;
    -tpr)
      shift #past argument
      TPR=$1
      shift #past value
      ;;
    -mdp)
      shift #past argument
      MDP=$1
      shift #past value
      ;;
    -start)
      shift #past argument
      RESTART=$1
      shift #past value
      ;;
    -o|--output)
      shift #past argument
      OUT=$1
      shift #past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
    *)
      echo "Unknown value $1"
      return 1
      ;;
  esac
done


echo 0 | gmx trjconv -f $TRAJ -o ${OUT}/AA-traj-${iso}.xtc \
  -s $TPR -n ../mapp-traj/index.ndx -pbc whole

python hack_topol.py $TOP ${OUT}/topol2m-${iso}_hacked.top

gmx grompp -p ${OUT}/topol2m-${iso}_hacked.top -f $MDP -c $RESTART -o ${OUT}/AA-COG-${iso}.tpr \
 -n ../structure_files/index-${iso}.ndx -maxwarn 1
rm mdout.mdp # clean-up

no_of_beads=$(grep "\[" ../structure_files/mapping.ndx | wc -l)
no_of_beads_minus_1=$( python -c "print( $no_of_beads - 1)" )

seq 0 $no_of_beads_minus_1 | gmx traj -f ${OUT}/AA-traj-${iso}.xtc -s ${OUT}/AA-COG-${iso}.tpr -oxt ${OUT}/mapped-${iso}.xtc \
  -n ../structure_files/mapping.ndx  -ng ${no_of_beads} -com


seq 0 $no_of_beads_minus_1 | gmx traj -f ${OUT}/AA-traj-${iso}.xtc -s ${OUT}/AA-COG-${iso}.tpr \
                                      -n ../structure_files/mapping.ndx -oxt ${OUT}/CG_motor2m-${iso}.gro -ng ${no_of_beads} -com -e 0

# Generate the CG tpr
gmx grompp -f mdp/CG_run.mdp -c ${OUT}/CG_motor2m-${iso}.gro -p ../structure_files/CG_system.top -o ${OUT}/CG-${iso}.tpr -maxwarn 1



rm \#*
