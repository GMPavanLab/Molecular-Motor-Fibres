#!/bin/bash


### Bonds and constraints 
traj=$1 #Location of trajectory
tpr=$2 #Location of tpr file
loc=$3 #Location of directory to store results
iso=R

DIR=${loc}/bonds_mapped
rm -rf $DIR/*
mkdir -p $DIR

NBONDS=$(grep "\[" ../structure_files/bonds.ndx | wc -l)
echo $NBONDS
IBOND=0
while [ $IBOND -lt $NBONDS ]; do
    echo $IBOND | gmx distance -f $traj -n ../structure_files/bonds.ndx -s $tpr -oall $DIR/bond_$IBOND.xvg -xvg none &> $DIR/temp_bonds.txt
    echo "---- bond"$IBOND" ----" >> $DIR/data_bonds.txt
    awk '/Average distance/ {print $3} /Standard deviation/ {print $3}' $DIR/temp_bonds.txt >> $DIR/data_bonds.txt
    #gmx analyze -f $DIR/bond_$IBOND.xvg -dist $DIR/distr_bond_$IBOND.xvg -xvg none -bw 0.001 
    rm -rf \#*
    let IBOND=$IBOND+1
done

# Clean-up
rm $DIR/temp_bonds.txt $DIR/\#*





### Angles

DIR=${loc}/angles_mapped
rm -rf $DIR/*
mkdir $DIR

NANGLES=$(grep "\[" ../structure_files/angles.ndx | wc -l)
echo $NANGLES
IANG=0
while [ $IANG -lt $NANGLES ]; do
    echo $IANG | gmx angle -f $traj -n ../structure_files/angles.ndx -ov $DIR/ang_$IANG.xvg -od $DIR/temp.xvg &> $DIR/temp_angles.txt
    echo "---- ang"$IANG" ----" >> $DIR/data_angles.txt
    awk '/< angle >/ {print $5} /Std. Dev./ {print $4}' $DIR/temp_angles.txt >> $DIR/data_angles.txt
    #gmx analyze -f $DIR/ang_$IANG.xvg -dist $DIR/distr_ang_$IANG.xvg -xvg none -bw 1.0 
    rm -rf \#*
    let IANG=$IANG+1
done

# Clean-up
rm $DIR/temp_angles.txt $DIR/temp.xvg $DIR/\#*

### Dihedrals

DIR=${loc}/dihedrals_mapped
rm -rf $DIR/*
mkdir $DIR

NDIHEDRALS=$(grep "\[" ../structure_files/dihedrals.ndx | wc -l)

IDIH=0
while [ $IDIH -lt $NDIHEDRALS ]; do
echo $IDIH
    echo $IDIH | gmx angle -type dihedral -f $traj -n ../structure_files/dihedrals.ndx -ov $DIR/dih_$IDIH.xvg -od $DIR/temp.xvg &> $DIR/temp_dihedrals.txt
    echo "---- dih"$IDIH" ----" >> $DIR/data_dihedrals.txt
    awk '/< angle >/ {print $5} /Std. Dev./ {print $4}' $DIR/temp_dihedrals.txt >> $DIR/data_dihedrals.txt
    #gmx analyze -f $DIR/dih_$IDIH.xvg -dist $DIR/distr_dih_$IDIH.xvg -xvg none -bw 1.0 
    rm -rf \#*
    let IDIH=$IDIH+1
done

# Clean-up
rm $DIR/temp_dihedrals.txt $DIR/temp.xvg $DIR/\#*

