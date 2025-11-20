#!/bin/bash

NCORES=16
GPU_ID=0
FILE_PATH=../structure_files/fibre_structure.gro
OUTPUT_NAME=tmp

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--file_name)
      shift #past argument
      FILE_PATH=$1
      shift
      ;;
    -p|--topology)
      shift #past argument
      TOPOLOGY_PATH=$1
      shift
      ;;
    -deffnm|--define_name)
       OUTPUT_NAME=$2
       shift;shift
       ;;
    -nt|--number_cores)
      NCORES=$2
      shift;shift
      ;;
    -gpu_id)
      GPU_ID=$2
      shift;shift
      ;;
    -nsteps|--number_steps)
      NSTEPS=$2
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
rm step*
rm \#*
rm ../structure_files/\#*
rm ../run/*/\#*
rm ../tmp/\#*

cp $TOPOLOGY_PATH ../tmp/topol.top
TOPOLOGY_PATH=../tmp/topol.top
cp $FILE_PATH ../tmp/box.gro

#echo 4 | gmx genrestr -f ../tmp/box.gro -o ../structure_files/posres_fibre.itp -fc 1000 1000 1000
printf "[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n" > ../structure_files/posres_fibre.itp
for i in {1..28}; do
  echo "$i 1 1000 1000 1000" >> ../structure_files/posres_fibre.itp
done

mkdir -p ../run/em_solvant
gmx grompp -p ${TOPOLOGY_PATH}  -c ../tmp/box.gro -f mdp/em_solvant.mdp \
  -r $FILE_PATH -o ../run/em_solvant/em.tpr  -maxwarn 1
gmx mdrun -deffnm ../run/em_solvant/em -gpu_id $GPU_ID -nt $NCORES

mkdir -p ../run/eq_solvant
gmx grompp -p ${TOPOLOGY_PATH} -c ../run/em_solvant/em.gro -f mdp/eq_solvant.mdp \
   -r $FILE_PATH  -o ../run/eq_solvant/eq.tpr  -maxwarn 1
gmx mdrun -deffnm ../run/eq_solvant/eq -gpu_id $GPU_ID -nt $NCORES

mkdir -p ../run/compress
gmx grompp -p ${TOPOLOGY_PATH} -c ../run/eq_solvant/eq.gro -f mdp/compress.mdp  -o ../run/compress/compress.tpr  -maxwarn 1
gmx mdrun -deffnm ../run/compress/compress -gpu_id $GPU_ID -nt $NCORES

#Expand box
gmx editconf -f ../run/compress/compress.gro -o ../tmp/distbox.gro -d 0.1
xlen=$(tail -n 1 ../tmp/distbox.gro | awk '{print $1}')
ylen=$(tail -n 1 ../run/compress/compress.gro | awk '{print $2}')
zlen=$(tail -n 1 ../tmp/distbox.gro | awk '{print $3}')
gmx editconf -f ../tmp/distbox.gro -o ../tmp/box_tight.gro -box $xlen $ylen $zlen 

mkdir -p ../run/compress2
gmx grompp -p ${TOPOLOGY_PATH} -c ../tmp/box_tight.gro -f mdp/compress2.mdp  -o ../run/compress2/compress.tpr  -maxwarn 1
gmx mdrun -deffnm ../run/compress2/compress -gpu_id $GPU_ID -nt $NCORES

cp ../run/compress2/compress.gro ../system/${OUTPUT_NAME}/compressed.gro
cp ../tmp/topol.top ../system/${OUTPUT_NAME}_topol.top