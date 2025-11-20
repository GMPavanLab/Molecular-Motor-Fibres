#!/bin/bash
#source /home/mattias/miniconda3/etc/profile.d/conda.sh
: ' Computes force field for a larger functional molecular molecule that needs to be split up into parts and then merged. The merging is performed using mainchain files. 
    Here, the subunits are unpolar residue (head), motor core + left leg[same leg as methyle group side] (core), and right leg (right_leg).
    Mainchain files has to be manually written, but for simillar functional motors, one can just ensure that dummy atoms on the subunits are named the same way as in the provided mainchain files.
  Inputs:
    -s|--sub_type: Subtype of motor. Used to find the correct mol2 file. The script matches the files core${SUBTYPE}-[RS].mol2 and right_leg${SUBTYPE}.mol2 in the molecular_muscle/structure_files/mol2 folder.
    -hc|--head_charge: Charge of the head.
    -cc|--core_charge: Charge of the core. 
    -lc|--leg_charge: Charge of the leg. 
  Returns: NONE
    Creates parametrised .gro file in molecular_muscle/structure_files/gro/motor${SUBTYPE}-[RS].gro
    Creates parametrised .pdb file in molecular_muscle/structure_files/pdb/motor${SUBTYPE}-[RS].pdb, oriented along the z-axis.
'
while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--sub_type)
      shift
      SUBTYPE=$1 
      shift
      ;;
    -hc|--head_charge)
      shift 
      HC=$1
      shift
      ;;
    -cc|--core_charge) 
      shift
      CC=$1
      shift
      ;;
    -lc|--leg_charge)
      shift
      LC=$1
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
     *)
      echo "Unknown option $1"
      return 1
  esac
done

for iso in R S; do

echo "Running ${SUBTYPE}${iso}-isomer"
mol2=../../structure_files/mol2
chain=../main_chain

PATH_CORE=${mol2}/core${SUBTYPE}B-${iso}.mol2 #Path residue. Note mol2 format
PATH_HEAD=${mol2}/head.mol2
PATH_RLEG=${mol2}/right_leg${SUBTYPE}.mol2

PATH_CORE_CHAIN=${chain}/core2m.chain #Path to mainchain files
PATH_HEAD_CHAIN=${chain}/head.chain
PATH_RLEG_CHAIN=${chain}/right_leg2m.chain

#Mol2_to_prepi generates .prepi and .frcmod files for each residue, with correct mainchain labeling.
echo "${PATH_HEAD} to prepi" 
. mol2_to_prepi.sh --file_name $PATH_HEAD --name HEAD --main_chain $PATH_HEAD_CHAIN --residue_name HED -nc $HC
echo "${PATH_RLEG} to prepi"
. mol2_to_prepi.sh --file_name $PATH_RLEG --name RLEG${SUBTYPE} --main_chain $PATH_RLEG_CHAIN --residue_name RLG -nc $LC
echo "${PATH_CORE} to prepi"
. mol2_to_prepi.sh --file_name $PATH_CORE --name CORE${SUBTYPE}-${iso} --main_chain $PATH_CORE_CHAIN --residue_name COR -nc $CC

: 'Build sequence joins the different residues and converts the tleap parameter outputfiles to .top and .gro files.
Order of residue names is not important. Note that they must correspond to a --name in mol2_to_prepi.
Order of sequence residues (-s) IS important. The first residue must ONLY have a tail_name defined, the last must ONLY have a head_name defined and any in the middle must have BOTH a head_name and a tail_name defined. (As well as potential mainchain atoms).
Results are sent to ../output.
'
echo "Patching Motor-${iso}"
. build_sequence.sh HEAD RLEG${SUBTYPE} CORE${SUBTYPE}-${iso} -s "HED COR RLG" -o motor${SUBTYPE}-${iso}

conda activate molecular-motors

printf "0\n" | gmx editconf -f ../output/motor${SUBTYPE}-${iso}.gro -o ../../structure_files/pdb/motor${SUBTYPE}-${iso}.pdb -box 10 10 4 -c -princ
cp ../output/motor${SUBTYPE}-${iso}.gro ../../structure_files/gro/motor${SUBTYPE}-${iso}.gro
cp ../output/motor${SUBTYPE}-${iso}.top  ../structure_files/top/motor${SUBTYPE}-${iso}.top
done

conda deactivate molecular-motors
