#!/bin/bash
#Makes prepi file from mol2 file and main_chain file. Used in combination with build_sequence command.

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--file_name) #Path to mol2 file
      shift #past argument
      FILE_NAME=$1
      shift
      ;;
    -n|--name) #What the files should be named
      shift #past argument
      NAME=$1
      shift
      ;;
    -nc|--net_charge) #What the files should be named
      shift #past argument
      NC=$1
      shift
      ;;
    -m|--main_chain) #Path to main_chain file
      shift #past argument
      CHAIN_NAME=$1
      shift
      ;;
    -rn|--residue_name) #What the prepi residue should be called
      shift #past argument
      RES_NAME=$1
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
  esac
done

mkdir -p ../logs/ #to store log-files
mkdir -p ../tmp/ #to store intermediate files
#Computes charges using bcc method. Results stored in .ac file.
antechamber -i $FILE_NAME -fi mol2 -o ../tmp/${NAME}.ac -fo ac -c bcc -at gaff2 -nc $NC -pf y 

#Generates the prepi file with correct mainchain labeling. Mainchain is defined by -m command. Note RES_NAME corresponds to the name later supplied to the sequence command.
prepgen -i ../tmp/${NAME}.ac -o ../tmp/${NAME}.prepi -f prepi -m ${CHAIN_NAME} -rn $RES_NAME -rf ${RES_NAME}.res

#Computes frcmod file with potential extra parameters. 
parmchk2 -i ../tmp/${NAME}.prepi -f prepi -o ../tmp/${NAME}.frcmod -s gaff2


