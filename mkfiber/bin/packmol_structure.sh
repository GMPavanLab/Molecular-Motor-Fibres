#!/bin/bash
#Generates a geometry with motors in cylinder. 4 motors in height, $N motors in circonference. 
#Requires packmol and gromacs

FILE_PATHS=()
MOLECULE_NUMBER=()
TOTAL_MOLECULES=0
NEUTRAL=0
mkdir -p ../tmp/
mkdir -p ../logs/

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--file_names)
      shift #past argument
      while [[ $1 == *.pdb ]]; do
         FILE_PATHS+=("$1")
         shift #past value
      done
      ;;
    -n|--number_of_molecules)
      shift #past argument
      while [[ $1 =~ ^[0-9]+$ ]]; do
         MOLECULE_NUMBER+=("$1")
	 TOTAL_MOLECULES=$((TOTAL_MOLECULES+$1))
         shift #past value
      done
      ;;
    -o|--output)
       OUTPUT_FILE_NAME=$2
        x=${2%.pdb}
        NAME=("${x##*/}")
       shift;shift
       ;;
    -neutral)
       NEUTRAL=1
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

CONSTRAIN_INNER=20 
CONSTRAIN_OUTER=$(($CONSTRAIN_INNER+20))
BOXLEN=$(($TOTAL_MOLECULES / 2))
echo "TOTAL NUMBER OF MOLECULES: $(($TOTAL_MOLECULES))"
echo "INNER RADIUS: $CONSTRAIN_INNER"
echo "OUTER RADIUS: $CONSTRAIN_OUTER"

cat > 'packmol.inp' << EOF

# All the atoms from diferent molecules will be at least 11 Angstrons
# apart at the solution
seed -1
tolerance 1.5
movebadrandom
# The output file name

output $OUTPUT_FILE_NAME
filetype pdb
writeout -1
EOF

for (( i=0; i<${#FILE_PATHS[@]}; i++));do
cat >> 'packmol.inp' << EOF

# Cylinder with ${MOLECULE_NUMBER[$i]} molecules in circomference and 4 molecules long.
structure ${FILE_PATHS[$i]}
  number $((${MOLECULE_NUMBER[$i]}))
  inside box -100 -100 -$BOXLEN 100 100 $BOXLEN

  atoms 21
    inside cylinder 0. 0. -200. 0. 0. 1. $CONSTRAIN_INNER. 400.
  end atoms

  atoms 72
    outside cylinder 0. 0. -200. 0. 0. 1. $CONSTRAIN_OUTER. 400.
  end atoms

  atoms 39 49 41 42 43 44 45 46 52 53 54 55 63 64 65 66 67 68 69 70 71 72 110 111 112 113 117 118
    radius 1.4

  constrain_rotation y 0. 0
  constrain_rotation x 0. 0
end structure 
EOF
done

if [ $NEUTRAL -eq 1 ];then
echo "Adding $((2*$TOTAL_MOLECULES)) Na+ to system"
cat >> 'packmol.inp' << EOF

# Adding $((2*$TOTAL_MOLECULES)) Na+ to system.
structure ../structure_files/NA.pdb
  number $((2*$TOTAL_MOLECULES))
  inside box -100 -100 -$BOXLEN 100 100 $BOXLEN

  atoms 1
    outside cylinder 0. 0. -200. 0. 0. 1. $(($CONSTRAIN_OUTER+5)). 400.
    inside cylinder 0. 0. -200. 0. 0. 1. $(($CONSTRAIN_OUTER+15)). 400.
  end atoms
end structure 
EOF
fi

#Run packmol on file
packmol < packmol.inp > ../logs/packmol_${NAME}.out
