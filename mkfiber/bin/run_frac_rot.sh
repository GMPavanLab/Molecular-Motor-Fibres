#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

: '
Utillity script that executes make_fiber.sh given a specific subtype of the MM, and the total number of initiall motors 
  - given an equal number of each isomer. Also computes center of polar chains for radius calculation.
'
GPU_ID=1
NCORES=16
N_SCALES=9

while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--subtype)
      shift #z
      SUBTYPE=$1
      shift
      ;;
    -gpu_id)
        shift #past argument
        GPU_ID=$1
        shift
        ;;
    -frac|--frac_rot)
        shift #past argument
        FRAC=$1
        shift
        ;;
    -n_scales)
        shift #past argument
        N_SCALES=$1
        shift
        ;;
    -deffnm)
        shift #past argument
        NAME=$1
        shift
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

val_B=$((50*$FRAC/100))
val_A=$((50-$val_B))
NATOMS=$(grep -c "ATOM" ../../structure_files/pdb/motor${SUBTYPE}-R.pdb)
TOPOLOGIES=()
PDBS=()
NRS=()
echo "Generating ${NAME}"

#If val_A > 0 include the A isomer
if [ $val_A -gt 0 ]; then
    TOPOLOGIES+=("../../structure_files/top/motor${SUBTYPE}-R.top")
    TOPOLOGIES+=("../../structure_files/top/motor${SUBTYPE}-S.top")
    PDBS+=("../../structure_files/pdb/motor${SUBTYPE}-R.pdb")
    PDBS+=("../../structure_files/pdb/motor${SUBTYPE}-S.pdb")
    NRS+=($val_A)
    NRS+=($val_A)
fi
#If val_B > 0 include the B isomer
if [ $val_B -gt 0 ]; then
    TOPOLOGIES+=("../../structure_files/top/motor${SUBTYPE}B-R.top")
    TOPOLOGIES+=("../../structure_files/top/motor${SUBTYPE}B-S.top")
    PDBS+=("../../structure_files/pdb/motor${SUBTYPE}B-R.pdb")
    PDBS+=("../../structure_files/pdb/motor${SUBTYPE}B-S.pdb")
    NRS+=($val_B)
    NRS+=($val_B)
fi
echo "PDBS: ${PDBS[@]}"
echo "TOPOLOGIES: ${TOPOLOGIES[@]}"
echo "NRS: ${NRS[@]}"
. make_fiber.sh -pdb ${PDBS[@]} -top ${TOPOLOGIES[@]} -n ${NRS[@]} -deffnm ${NAME} -n_scales $N_SCALES -gpu_id $GPU_ID
  
