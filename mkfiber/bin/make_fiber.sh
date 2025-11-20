#!/bin/bash
: '
MASTER RUNFILE THAT GENERATES ANY FIBER STRUCTURE
Genereates one fibre from a list of motor structures and number of each respective molecule. Greatly limited by the initial packmol packageing. 
If your system is unstable, try running with a different seed. If it is still unstable, try running with a different number of motors.
Also try ensuring that the motors are as "straight" as possible to facilitate packing. 

INPUTS:
  -pdb: paths
    list of pdb files of motors. NOTE: Oriented to principal axis
  -top: paths
    list of topology files of respective motors.
  -n: int > 0
    list of number of each molecule. NOTE: Sum between 80 and 110 for best performance. 
  -deffnm: string
    name of the output files.
  -gpu_id: int
    gpu id to use for gromacs.

OUTPUTS:
  -../system/${DEFFNM}/: folder with all the output files.
      .xtc, .cpt, .gro and .edr files of the final NPT.
      .pdb structure without water, and .top file of the final solved fibre.
  -../output/NPT/(nvt|npt|npt_long)_${DEFFNNM}_(z|density|pressure|temperature).xvg: xvg files form the many equillibration steps.

TO DO:
  In the current version, the program assumes each motor has the charge -2, but it would not be very difficult to simply 
    change the number of ions added in make_fibre_structure.sh if that would be needed.
  Would be nice to increase length of equillibration to 50 - 100 bs, but this would require a lot of computational time.

EX:
  . make_fiber.sh -pdb ../../structure_files/pdb/motor2m-R.pdb ../../structure_files/pdb/motor2m-S.pdb\
    -top ../../structure_files/top/motor2m.top ../../structure_files/top/motor2m-S.top\
    --number_of_molecules 55 55 -deffnm example -gpu_id 0
'


NCORES=16
PDBS=()
TOPS=()
NRS=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -n|--number_of_molecules)
      shift 
      while [[ $1 =~ ^[0-9]+$ ]]; do
         NRS+=("$1")
         shift
      done
      ;;
    -pdb|--pdb_files)
      shift 
      while [[ $1 == *.pdb ]]; do
         PDBS+=("$1")
         shift 
      done
      ;;
    -top |--top_files)
      shift 
      while [[ $1 == *.top ]]; do
         TOPS+=("$1")
         shift
      done
      ;;
    -deffnm)
      DEFFNM=$2
      shift;shift
      ;;
    -n_scales)
      N_SCALES=$2
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

mkdir -p ../tmp/
mkdir -p ../output/${OUTPUT_PATH}
mkdir -p ../system/${DEFFNM}
#Packmol to generate initial structure
. packmol_structure.sh -f ${PDBS[@]} -n ${NRS[@]} -o ../system/${DEFFNM}/${DEFFNM}.pdb -neutral
#Compress fiber through energy minimization
TOPOLOGY_PATH=../tmp/${DEFFNM}.top
bash compress_fibre.sh -f ../system/${DEFFNM}/${DEFFNM}.pdb -p $TOPOLOGY_PATH\
    -t  ${TOPS[@]}  -n ${NRS[@]} -gpu_id $GPU_ID -n_scales $N_SCALES

#Duplicate size of fiber, solvate in water, and equilibriate
bash equilibriate_large_fiber.sh -f ../output/${DEFFNM}/${DEFFNM}_dense.gro -p $TOPOLOGY_PATH \
    -n 2 -deffnm $DEFFNM -nt $NCORES -gpu_id $GPU_ID
