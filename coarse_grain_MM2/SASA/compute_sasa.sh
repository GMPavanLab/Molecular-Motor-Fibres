#!/bin/bash
: '
Measures SASA of a molecule in a trajectory.
INPUT:
    $1: (MM2|MM2b-R)
        Molecule name
    $2: str
        version of the coarse-grained model. Must have a trajectory file in ../match-dist/${version}/run/. Not needed if $3 is "AA", then just write any str.
    $3: (AA|CG)
        Atomistic or Coarse-Grained
OUTPUT:
    (AA|CG)/(MM2|MM2b-R)_SASA.xvg
        SASA of the molecule in the trajectory
'
molecule=$1
version=$2
on=$3
size_1=${#molecule}

if [ $on == "AA" ]; then
# Atomistic
mkdir -p AA
cd       AA
cp ../vdwradii_AA.dat  vdwradii.dat
# Compute SASA
echo 1 | gmx sasa -f ../../mapp-traj/AA-traj-${molecule}.xtc -s ../../mapp-traj/AA-COG-${molecule}.tpr -ndots 4800 -probe 0.185  -o ${molecule}_SASA-AA.xvg
#gmx sasa -f ../../mapp-traj/AA-traj-${molecule}.xtc -s ../../mapp-traj/AA-COG-${molecule}.tpr -n ../sasa_index_MM2b-R.ndx -ndots 4800 -probe 0.185 -o ${molecule}_SASA_AA-STATOR.xvg

gmx analyze -f ${molecule}_SASA-AA.xvg -bw 0.01 -dist ${molecule}_distr-SASA-AA.xvg &> ${molecule}_temp_sasa.txt
awk '/SS1/ {print "average: " $2}/SS1 / {print "st.dev.: " $3}' ${molecule}_temp_sasa.txt > ${molecule}_data_sasa_AA.xvg

# Compute Connolly Surface
echo 1 | gmx trjconv -f ../../AA-reference/em/em_${molecule}.gro -o ${molecule}-AA-min.gro -s ../../AA-reference/em/em_${molecule}.tpr -pbc whole
echo 1 | gmx sasa -s ${molecule}-AA-min.gro -o temp.xvg -probe 0.185 -ndots 240 -q ${molecule}_surf-AA.pdb
rm \#* ${molecule}_temp_sasa.txt average.xvg errest.xvg temp.xvg
cd       ..
else

# Coarse-Grained
mkdir -p CG
cd       CG
cp ../vdwradii_CG.dat  vdwradii.dat
# Compute SASA
echo 2 | gmx sasa -f ../../match-dist/${version}/run/run_${molecule}.xtc -s ../../match-dist/${version}/run/run_${molecule}.tpr -ndots 4800 -probe 0.185  -o ${molecule}_SASA-CG.xvg
gmx analyze -f ${molecule}_SASA-CG.xvg -bw 0.01 -dist ${molecule}_distr-SASA-CG.xvg &> ${molecule}_temp_sasa.txt
awk '/SS1/ {print "average: " $2}/SS1 / {print "st.dev.: " $3}' ${molecule}_temp_sasa.txt > ${molecule}_data_sasa_CG.xvg

# Compute Connolly Surface (first need to map the energy-minimized AA snapshot!)
no_of_beads=$(grep "\[" ../../structure_files/mapping_${molecule}.ndx | wc -l)
no_of_beads_minus_1=$( python -c "print( $no_of_beads - 1)" )

seq 0 $no_of_beads_minus_1 | gmx traj -f ../AA/${molecule}-AA-min.gro -s ../../mapp-traj/AA-COG-${molecule}.tpr \
                                      -n ../../structure_files/mapping_${molecule}.ndx -oxt ${molecule}-MAPPED-min.gro -ng ${no_of_beads} -com 
echo 2 | gmx trjconv -f ${molecule}-MAPPED-min.gro -o ${molecule}-MAPPED-min-CGatomnames.gro -s ../../match-dist/${version}/eq5/eq_${molecule}.tpr
echo 2 | gmx sasa -s ${molecule}-MAPPED-min-CGatomnames.gro -o temp.xvg -probe 0.185 -ndots 240 -q ${molecule}_surf-CG.pdb
rm \#* ${molecule}_temp_sasa.txt average.xvg errest.xvg temp.xvg
cd       ..

fi