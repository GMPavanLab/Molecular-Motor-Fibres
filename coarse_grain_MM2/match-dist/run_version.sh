#!/bin/bash
: '
Runs the CG model for the specified molecule and topology version. Measures the parameters provided in ../structure_files/(bonds|angles|dihedrals).ndx.
INPUTS:
    -v|--version: dir
        version of the topology to use. Must be a folder in the current directory, containing a topology file CG_system_base.top
    -mol|--molecule: (MM2 | MM2b-R)
        molecule to simulate (MM2 or MM2b-R)
    -nt|--number_cores: number of cores to use
    -gpu_id: GPU ID to use
OUTPUTS:
    -${version}/(angles|bonds|dihedrals)_mapped/* :
        directories containing the measured parameters as xvg files for distributions and .txt files with computed averages and standard deviations.
'

NCORES=16
GPU_ID=0

while [[ $# -gt 0 ]]; do
  case $1 in
    -v|--version)
      shift #past argument
      version=$1
      shift #past value
      ;;
    -mol|--molecule)
      shift #past argument
      MOL=$1 #MM2 or MM2b-R
      shift #past value
      ;;
    -nt|--number_cores)
       NCORES=$2
       shift;shift
       ;;
    -gpu_id)
        GPU_ID=$2
        shift;shift
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

solute=../mapp-traj/CG_${MOL}.gro
solvent_box=../structure_files/box_CG_W_eq.gro
solvent_name=W
solvent_atoms=1

cp ${version}/CG_system_base.top ${version}/CG_system.top
gmx solvate -cp ${solute} -cs ${solvent_box} -p ${version}/CG_system.top -o ${version}/initial_${MOL}.gro -box 5 5 5

echo "Adding ions"
gmx grompp -f mdp/ions.mdp -c ${version}/initial_${MOL}.gro -p ${version}/CG_system.top -o ${version}/ions.tpr -maxwarn 1
printf "6\n" | gmx genion -s ${version}/ions.tpr -o ${version}/initial_${MOL}.gro -p ${version}/CG_system.top -pname NA -neutral

mkdir -p ${version}/em
gmx grompp -p ${version}/CG_system.top -c ${version}/initial_${MOL}.gro -f mdp/em.mdp  -o ${version}/em/em_${MOL}.tpr -po ${version}/em/em.mdp  -maxwarn 1
gmx mdrun -v -deffnm ${version}/em/em_${MOL} -gpu_id $GPU_ID -nt $NCORES 

mkdir -p ${version}/output
printf "8\n0\n" | gmx energy -f ${version}/em/em_${MOL}.edr -o ${version}/output/em_${MOL}.xvg

mkdir -p ${version}/eq1
gmx grompp -p ${version}/CG_system.top -c ${version}/em/em_${MOL}.gro -f mdp/eq1.mdp  -o ${version}/eq1/eq_${MOL}.tpr -po ${version}/eq1/eq.mdp  -maxwarn 1
gmx mdrun -deffnm ${version}/eq1/eq_${MOL} -gpu_id $GPU_ID -nt $NCORES -bonded gpu

for i in {2..5}; do
    mkdir -p ${version}/eq$i
    gmx grompp -p ${version}/CG_system.top -c ${version}/eq$(( $i - 1 ))/eq_${MOL}.gro -r ${version}/eq$(( $i - 1 ))/eq_${MOL}.gro -f mdp/eq${i}.mdp \
     -o ${version}/eq${i}/eq_${MOL}.tpr -po ${version}/eq${i}/eq.mdp  -maxwarn 1
    gmx mdrun -v -deffnm ${version}/eq${i}/eq_${MOL} -gpu_id $GPU_ID -nt $NCORES  -bonded gpu
    printf "12\n0\n" | gmx energy -f ${version}/eq${i}/eq_${MOL}.edr -o ${version}/output/eq${i}-temp_${MOL}.xvg
    printf "13\n0\n" | gmx energy -f ${version}/eq${i}/eq_${MOL}.edr -o ${version}/output/eq${i}-pressure_${MOL}.xvg
    printf "19\n0\n" | gmx energy -f ${version}/eq${i}/eq_${MOL}.edr -o ${version}/output/eq${i}-density_${MOL}.xvg
done

mkdir -p ${version}/run

gmx grompp -p ${version}/CG_system.top -c ${version}/eq5/eq_${MOL}.gro -f mdp/run.mdp \
 -o ${version}/run/run_${MOL}.tpr -po ${version}/run/run_${MOL}.mdp 
gmx mdrun -v -deffnm ${version}/run/run_${MOL} -gpu_id $GPU_ID -nt $NCORES -bonded gpu

rm ${version}/*/\#*
. compute_distr.sh ${version}/run/run_${MOL}.xtc ${version}/run/run_${MOL}.tpr $version