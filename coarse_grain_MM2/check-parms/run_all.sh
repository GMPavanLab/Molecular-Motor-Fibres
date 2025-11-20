#!/bin/bash
: '
Measures parameters in the AA reference for comparison with literature values. The results are then analysed using analysis.ipynb
OUTPUT:
    -output/* : .xvg - files for 4 critical angles, 1 bond, and 2 dihedrals.
'
mkdir -p referece
cp ../AA-reference/run/md_MM2.gro referece/MM2-AA.gro
traj=../AA-reference/run/md_MM2.xtc
tpr=../AA-reference/run/md_MM2_run.tpr
DIR=output
rm -r output
mkdir -p output

NBONDS=$(grep "\[" bonds-AA.ndx | wc -l)

for ((IBOND=0; IBOND<$NBONDS; IBOND++)); do
    echo "BOND "$IBOND
    echo $IBOND | gmx distance -f $traj -n bonds-AA.ndx -s $tpr -oall $DIR/bond_$IBOND.xvg -xvg none &> $DIR/temp_bonds.txt
    echo "---- bond"$IBOND" ----" >> $DIR/data_bonds.txt
    awk '/Average distance/ {print $3} /Standard deviation/ {print $3}' $DIR/temp_bonds.txt >> $DIR/data_bonds.txt
    rm -rf \#*
done

rm $DIR/temp_bonds.txt $DIR/temp.xvg $DIR/\#*

NANGLES=$(grep "\[" angles-AA.ndx | wc -l)
for ((IANG=0; IANG<$NANGLES; IANG++)); do
    echo "ANGLE "$IANG
    echo $IANG | gmx angle -f $traj -n angles-AA.ndx -ov $DIR/ang_$IANG.xvg -od $DIR/temp.xvg &> $DIR/temp_angles.txt
    echo "---- ang"$IANG" ----" >> $DIR/data_angles.txt
    awk '/< angle >/ {print $5} /Std. Dev./ {print $4}' $DIR/temp_angles.txt >> $DIR/data_angles.txt
    rm -rf \#
done

# Clean-up
rm $DIR/temp_angles.txt $DIR/temp.xvg $DIR/\#*





### Dihedrals

NDIHEDRALS=$(grep "\[" dihedrals-AA.ndx | wc -l)

for ((IDIH=0; IDIH<$NDIHEDRALS; IDIH++)); do
    echo "DIHEDRAL "$IDIH
    echo $IDIH | gmx angle -type dihedral -f $traj -n dihedrals-AA.ndx -ov $DIR/dih_$IDIH.xvg -od $DIR/temp.xvg &> $DIR/temp_dihedrals.txt
    echo "---- dih"$IDIH" ----" >> $DIR/data_dihedrals.txt
    awk '/< angle >/ {print $5} /Std. Dev./ {print $4}' $DIR/temp_dihedrals.txt >> $DIR/data_dihedrals.txt
    #gmx analyze -f $DIR/dih_$IDIH.xvg -dist $DIR/distr_dih_$IDIH.xvg -xvg none -bw 1.0 
    rm -rf \#*
done

# Clean-up
rm $DIR/temp_dihedrals.txt $DIR/temp.xvg $DIR/\#*

