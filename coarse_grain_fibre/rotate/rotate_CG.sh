#!/bin/bash
#WORKS USING RANDOM SEED  -1512343793
NCORES=$1
GPU_ID=$2
KAPPA=$3
NROT=$4
NAME=$5
MOL=MM2
NMOLS=220
NATOMS=28
STEPLEN=5000
START=1000

cp ../CG-reference/eq5/eq.gro FMM.gro
cp ../CG-reference/eq5/eq.cpt FMM.cpt
cp ../CG-reference/index.ndx index.ndx
cp ../CG-reference/CG_system.top topol.top
sed -i '/; Include chain topologies/a #include "../structure_files/CG_motorb-R.itp"' topol.top
sed -i '/; Include chain topologies/a #include "../structure_files/CG_motorb-S.itp"' topol.top
sed -i '/MOTOR-R/ i\MOTORB-R '"$NROT"'' topol.top
sed -i '/MOTOR-S/ i\MOTORB-S '"$NROT"'' topol.top
sed -i -E '/MOTOR-[RS]/s/[0-9]+/'"$((55 - $NROT))"'/' topol.top

rm step*
rm run/\#*
mkdir -p COLVARS
rm -f COLVARS/${NAME}*
#Measure state of one molecule MD
cat > plumed.dat << EOF
dist_A: RMSD REFERENCE=CG_MM2a_lit.pdb TYPE=OPTIMAL
dist_B: RMSD REFERENCE=CG_MM2b_lit.pdb TYPE=OPTIMAL 
dist_I1: RMSD REFERENCE=CG_MM2I1_lit.pdb TYPE=OPTIMAL

EOF

#Steered MD
for j in 0 55 110 165; do
    for ((i=0; i<$(( $NROT )); i+=1)); do
    factor=$(( $i * $NATOMS + $j*$NATOMS ))
    RSp=""
    RSm="-"
    RSL="L"
    #If j == 55 or 165, RSp = "-" RSm = ""
    if [ $j -eq 55 ] || [ $j -eq 165 ]; then
        RSp="-"
        RSm=""
        RSL="U"
    fi

        cat >> 'plumed.dat' << EOF
alpha${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor + 7 )) AXIS=$(( $factor + 1 )),$(( $factor + 3 )) VECTOR2=$(( $factor + 8 )),$(( $factor + 14 ))
beta${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 1 )),$(( $factor + 14 ))
alpha2_${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +7 )) AXIS=$(( $factor + 1 )),$(( $factor +3 )) VECTOR2=$(( $factor + 8 )),$(( $factor + 13 ))
beta2_${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 1 )),$(( $factor + 13 ))
gamma${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 1 )),$(( $factor + 8 ))
gamma2_${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 13 )),$(( $factor + 14 ))


#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint${i}_${j}: MOVINGRESTRAINT ...
    ARG=alpha${i}_${j},alpha2_${i}_${j},beta${i}_${j},beta2_${i}_${j},gamma${i}_${j},gamma2_${i}_${j} VERSE=B,B,B,${RSL},${RSL},B
    AT0=${RSp}0.5,${RSp}0.5,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP0=$START                   KAPPA0=0,0,0,0,0,0
    AT1=${RSp}0.5,${RSp}0.5,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP1=$(($START + 1000))       KAPPA1=$KAPPA,$KAPPA,$KAPPA,$KAPPA,40,50
    AT2=${RSp}2.8,${RSp}2.8,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP2=$(($START + 1000 + 2*$STEPLEN))       KAPPA2=$KAPPA,$KAPPA,$KAPPA,$KAPPA,40,50
...
EOF
done
done

#Compute radius
for ((i=0; i<$NMOLS; i+=1)); do
    cat >> 'plumed.dat' << EOF

hed$i: CENTER ATOMS=$(( $i*$NATOMS + 25 ))-$(( $i*$NATOMS + 28 ))
core$i: CENTER ATOMS=$(( $i*$NATOMS + 1 ))-$(( $i*$NATOMS + 14 ))
tail$i: CENTER ATOMS=$(( $i*$NATOMS + 15 ))-$(( $i*$NATOMS + 24 ))
oxy$i: CENTER ATOMS=$(( 17 + $i * $NATOMS )),$(( 18 + $i * $NATOMS ))
EOF
done

#Dump positions for radii calculation
for name in hed core tail oxy; do
    printf "DUMPATOMS STRIDE=100 FILE=COLVARS/${NAME}_${name}.xyz ATOMS=${name}0" >> plumed.dat
    for ((i=1; i<$(($NMOLS)); i+=1)); do
        printf ",${name}$i" >> plumed.dat
    done
    printf "\n" >> plumed.dat
done

#Save states for first molecule
printf "PRINT STRIDE=100 ARG=dist_A,dist_I1,dist_B FILE=COLVARS/${NAME}_COLVAR\n" >> plumed.dat

#Run MD
mkdir -p run
gmx grompp -p topol.top -n index.ndx -f mdp/run-CG.mdp -c FMM.gro -t FMM.cpt -o run/run.tpr -maxwarn 1
gmx mdrun -nt $NCORES -gpu_id $GPU_ID -deffnm run/run -plumed plumed.dat -v
