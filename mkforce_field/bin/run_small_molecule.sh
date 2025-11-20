#!/bin/bash
# Generates a force field for smaller molecule that does not require being split into parts and then merged. 
# Here only used to compute force field of the motor cores used for coarse graining.

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
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
     *)
      echo "Unknown option $1"
      return 1
  esac
done

#conda activate amber

antechamber -i $FILE_NAME.mol2 -fi mol2 -o ../tmp/${NAME}.mol2 -fo mol2 -c bcc -at gaff2 -nc 0
parmchk2 -i ../tmp/${NAME}.mol2 -f mol2 -o ../tmp/${NAME}.frcmod -s gaff2

#Load parameters
cat > "tleap.in" << EOF
source leaprc.gaff2
MOL = loadmol2 ../tmp/${NAME}.mol2
loadamberparams ../tmp/${NAME}.frcmod 
saveamberparm MOL ../tmp/${NAME}.prmtop ../tmp/${NAME}.inpcrd
quit
EOF

tleap -s -f tleap.in > ../logs/tleap.out #Run tleap

python - << END
import parmed as pmd
parm = pmd.load_file('../tmp/${NAME}.prmtop', '../tmp/${NAME}.inpcrd')
parm.save('../output/${NAME}.top', format='gromacs')
parm.save('../output/${NAME}.gro')
END

