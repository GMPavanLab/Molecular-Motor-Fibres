#!/bin/bash

molecule=MOTOR
version=$1
size_1=${#molecule}
res=MOTOR

if [[ $version == "AA" ]]; then
    # Atomistic
    mkdir -p AA
    cd       AA
    cp ../vdwradii_AA.dat  vdwradii.dat
    # Compute SASA
    echo ${molecule} | gmx sasa -f ../../mapp-traj/AA-traj-R.xtc -n ../../structure_files/index-R.ndx -s ../../mapp-traj/AA-COG-R.tpr -ndots 4800 -probe 0.185  -o SASA-AA.xvg
    gmx analyze -f SASA-AA.xvg -bw 0.01 -dist distr-SASA-AA.xvg &> temp_sasa.txt
    awk '/SS1/ {print "average: " $2}/SS1 / {print "st.dev.: " $3}' temp_sasa.txt > data_sasa_AA.xvg

    # Compute Connolly Surface
    echo ${molecule} | gmx trjconv -f ../../AA-reference/em/em2m-R.gro -n ../../structure_files/index-R.ndx -o ${molecule}-AA-min.gro -s ../../AA-reference/em/em2m-R.tpr -pbc whole
    echo ${molecule} | gmx sasa -s ${molecule}-AA-min.gro -n ../../structure_files/index-R.ndx -o temp.xvg -probe 0.185 -ndots 240 -q surf-AA.pdb
    rm \#* temp_sasa.txt average.xvg errest.xvg temp.xvg
    cd       ..
fi

# Coarse-Grained
mkdir -p CG
cd       CG
cp ../vdwradii_CG.dat  vdwradii.dat
# Compute SASA
echo ${res} | gmx sasa -f ../../match-dist/${version}/run/run.xtc -n ../../match-dist/${version}/index.ndx -s ../../match-dist/${version}/run/run.tpr -ndots 4800 -probe 0.185  -o SASA-CG.xvg
gmx analyze -f SASA-CG.xvg -bw 0.01 -dist distr-SASA-CG.xvg &> temp_sasa.txt
awk '/SS1/ {print "average: " $2}/SS1 / {print "st.dev.: " $3}' temp_sasa.txt > data_sasa_CG.xvg

# Compute Connolly Surface (first need to map the energy-minimized AA snapshot!)
no_of_beads=$(grep "\[" ../../structure_files/mapping.ndx | wc -l)
no_of_beads_minus_1=$( python -c "print( $no_of_beads - 1)" )

seq 0 $no_of_beads_minus_1 | gmx traj -f ../AA/${molecule}-AA-min.gro -s ../../mapp-traj/AA-COG-R.tpr \
                                      -n ../../structure_files/mapping.ndx -oxt ${molecule}-MAPPED-min.gro -ng ${no_of_beads} -com 
echo ${res} | gmx trjconv -f ${molecule}-MAPPED-min.gro -n ../../match-dist/${version}/index.ndx -o ${molecule}-MAPPED-min-CGatomnames.gro -s ../../match-dist/${version}/eq5/eq.tpr
echo ${res} | gmx sasa -s ${molecule}-MAPPED-min-CGatomnames.gro -n ../../match-dist/${version}/index.ndx -o temp.xvg -probe 0.185 -ndots 240 -q surf-CG.pdb
rm \#* temp_sasa.txt average.xvg errest.xvg temp.xvg
cd       ..

