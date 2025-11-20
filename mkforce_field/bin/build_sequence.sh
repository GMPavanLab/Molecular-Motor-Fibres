#!/bin/bash
#Builds molecule from sequence of prepi-files. 

FILE_NAMES=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--sequence) #Sequence of prepi residues in format: "RS1 RS2 RS2"
      shift #past argument
      SEQUENCE=$1
      shift
      ;;
    -o|--output_file_name)
      shift #past argument
      OUTPUT_NAME=$1
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
    *)
      FILE_NAMES+=($1) #Adds name to different residue files
      shift
      ;;
  esac
done

mkdir -p ../output/ #To store results

#Generate tleap.in file
echo "source leaprc.gaff2" > tleap.in
for NAME in ${FILE_NAMES[@]}; do
#Load parameters
cat >> "tleap.in" << EOF
loadamberparams ../tmp/${NAME}.frcmod 
loadamberprep ../tmp/${NAME}.prepi
EOF
done
#Generate sequence and store resulting structure in pdb-file
cat >> "tleap.in" << EOF
mol = sequence {$SEQUENCE}
savemol2 mol ../output/${OUTPUT_NAME}.mol2 0
saveamberparm mol ../tmp/${OUTPUT_NAME}.prmtop ../tmp/${OUTPUT_NAME}.inpcrd
quit
EOF

tleap -s -f tleap.in > ../logs/tleap.out #Run tleap

rm ../output/${OUTPUT_NAME}.top
rm ../output/${OUTPUT_NAME}.gro

#Convert amber files to gromacs files and store them in output. 
python - << END
import parmed as pmd
parm = pmd.load_file('../tmp/${OUTPUT_NAME}.prmtop', '../tmp/${OUTPUT_NAME}.inpcrd')
parm.save('../output/${OUTPUT_NAME}.top', format='gromacs')
parm.save('../output/${OUTPUT_NAME}.gro')
END

