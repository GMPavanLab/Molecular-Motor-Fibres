#!/bin/bash
GPU_ID=0
NCORES=16

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
    -cion|--current_ion)
       ION_IN=$2
       shift;shift
       ;;
    -nion|--new_ion)
         ION_OUT=$2
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

cp $TOPOLOGY_PATH ../tmp/$OUTPUT_NAME.top
cp $FILE_PATH ../tmp/$OUTPUT_NAME.gro
TOPOLOGY_PATH=../tmp/$OUTPUT_NAME.top

#Just edit the ion type inline and pray the system doesnt explode
sed -i "s/$ION_IN/$ION_OUT/g" $TOPOLOGY_PATH
sed -i "s/$ION_IN/$ION_OUT/g" ../tmp/$OUTPUT_NAME.gro



mkdir -p ../run/em
gmx grompp -p ${TOPOLOGY_PATH} -c ../tmp/$OUTPUT_NAME.gro -f mdp/em.mdp  -o ../run/em/em.tpr  -maxwarn 1
gmx mdrun -v -deffnm ../run/em/em -gpu_id $GPU_ID -nt $NCORES

mkdir -p ../run/eq5
gmx grompp -p ${TOPOLOGY_PATH} -c ../run/em/em.gro -f mdp/eq5.mdp \
    -n $INDEX -o ../run/eq5/eq.tpr  -maxwarn 1
gmx mdrun -v -deffnm ../run/eq5/eq -gpu_id $GPU_ID -nt $NCORES
