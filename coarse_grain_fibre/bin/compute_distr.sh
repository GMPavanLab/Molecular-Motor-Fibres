#!/bin/bash


### Bonds and constraints 
traj=$1 #Location of trajectory
tpr=$2 #Location of tpr file
loc=$3 #Location of directory to store results
dn=$4 #Time between frames in ps
CG=$5 #1 if cg 0 if AA
NMOLS=200
NATOMS=28


mkdir -p ../colvars
cd ../colvars
[ "$CG" -eq 0 ] && echo 0 || echo 2 | gmx trjconv -f $traj -s $tpr -skip $dn -o subset.xtc
traj=subset.xtc

NBONDS=$(grep "\[" ../structure_files/bonds.ndx | wc -l)
python ../bin/expand_mapping.py ../structure_files/bonds.ndx bonds_expanded.ndx $NATOMS $NMOLS
IBOND=0
for (( i=0; i<$NMOLS; i++)); do
    echo "Bonds molecule $i"
    DIR=${loc}/MOTOR$i/bonds_mapped
    rm -rf $DIR/*
    mkdir -p $DIR
    while [ $IBOND -lt $(( $NBONDS * $i + $NBONDS )) ]; do
        echo $IBOND | gmx distance -f $traj -n bonds_expanded.ndx -s $tpr -oall $DIR/bond_$(( $IBOND - $i*$NBONDS )).xvg -xvg none &>> tmp.out
        rm -rf \#*
        rm -rf $DIR/\#*
        let IBOND=$IBOND+1
    done
done

### Angles
NANGLES=$(grep "\[" ../structure_files/angles.ndx | wc -l)
python ../bin/expand_mapping.py ../structure_files/angles.ndx angles_expanded.ndx $NATOMS $NMOLS
IANG=0
for (( i=0; i<$NMOLS; i++ )); do
    echo "Angles molecule $i"
    DIR=${loc}/MOTOR$i/angles_mapped
    rm -rf $DIR/*
    mkdir -p $DIR
    while [ $IANG -lt $(( $NANGLES * $i + $NANGLES )) ]; do
        echo $IANG | gmx angle -f $traj -n angles_expanded.ndx -ov $DIR/ang_$(( $IANG - $i*$NANGLES )).xvg -od $DIR/temp.xvg &>> tmp.out
        rm -rf \#*
        rm -rf $DIR/\#*
        let IANG=$IANG+1
    done    
done
### Dihedrals

NDIHEDRALS=$(grep "\[" ../structure_files/dihedrals.ndx | wc -l)
python ../bin/expand_mapping.py ../structure_files/dihedrals.ndx dihedrals_expanded.ndx $NATOMS $NMOLS
IDIH=0
for (( i=0; i<$NMOLS; i++ )); do
    echo "Dihedrals molecule $i"
    DIR=${loc}/MOTOR$i/dihedrals_mapped
    rm -rf $DIR/*
    mkdir -p $DIR
    while [ $IDIH -lt $(( $NDIHEDRALS *$i + $NDIHEDRALS )) ]; do
        echo $IDIH | gmx angle -type dihedral -f $traj -n dihedrals_expanded.ndx -ov $DIR/dih_$(( $IDIH - $i*$NDIHEDRALS )).xvg -od $DIR/temp.xvg &>> tmp.out
        rm -rf \#*
        rm -rf $DIR/\#*
        let IDIH=$IDIH+1
    done    
done
cd ../bin   