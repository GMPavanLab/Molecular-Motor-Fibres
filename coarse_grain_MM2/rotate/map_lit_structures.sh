#!/bin/bash/
#Mapps the structures from the literature for the different configurations to the CG structures
for structure in a b I1; do 
cat > ../structure_files/MM2${structure}_tmp2.gro << EOF
GROningen MAchine for Chemical Simulation
   49
EOF
    awk '{ printf "%8s %6s %4d %8.3f %8.3f %8.3f\n", "1OR2", $1, NR, -$2/10, $3/10, $4/10 }' ../structure_files/MM2${structure}_LIT.txt >> ../structure_files/MM2${structure}_tmp2.gro
    echo "   1.34130   1.54560   1.46670\n" >> ../structure_files/MM2${structure}_tmp2.gro
    echo "0" | gmx editconf -f ../structure_files/MM2${structure}_tmp2.gro -o ../structure_files/MM2${structure}_lit.pdb -princ
    seq 0 13 | gmx traj -f ../structure_files/MM2${structure}_lit.pdb -s ../structure_files/MM2${structure}_lit.pdb -n ../structure_files/mapping_lit${structure}.ndx -oxt ../structure_files/CG_MM2${structure}_lit.pdb -com -ng 14
    rm ../structure_files/MM2${structure}_tmp2.gro
    rm ../structure_files/\#*
done