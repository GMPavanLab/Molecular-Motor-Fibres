#!/bin/bash

mkdir -p ../run
cp ../structure_files/fibre.top ../run/topol.top
rm ../COLVARS/bck* ../positions/bck* step*
rm plumed.dat
rm ../run/\#*
NCORES=$1
GPU_ID=$2
KAPPA=$3
tropic=$4
NROT=$5
NMOLS=100
NATOMS=152
cat > 'plumed.dat' << EOF
MOLINFO STRUCTURE=../structure_files/fibre_structure.pdb
EOF

printf 'WHOLEMOLECULES ENTITY0=1' >> plumed.dat

#Make motors whole 
for ((i=2; i<=$(($NMOLS*$NATOMS)); i+=1)); do
    printf ",$i" >> plumed.dat
done
for ((i=$(( $NMOLS*$NATOMS + 2*$NMOLS +1 )); i<=$(( $NMOLS*$NATOMS*2 + 2*$NMOLS )); i+=1)); do
    printf ",$i" >> plumed.dat
done
printf "\n" >> plumed.dat

#Perform steered metad on dihedrals
j=0
    for ((i=0; i<$(( $NROT )); i+=1)); do
    factor=$(( $i * $NATOMS + $j ))
        cat >> 'plumed.dat' << EOF
com5_${i}_${j}: COM ATOMS=$((117+$factor)),$((118+$factor))
com1_${i}_${j}: COM ATOMS=$((66+$factor)),$((67+$factor)),$((109+$factor))
com8_${i}_${j}: COM ATOMS=$((63+$factor)),$((64+$factor))
com13_${i}_${j}: COM ATOMS=$((43+$factor)),$((44+$factor)),$((45+$factor)),$((50+$factor))
phi${i}_${j}: TORSION ATOMS=$(( 63 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 110 + $j ))
psi${i}_${j}: TORSION ATOMS=$(( 55 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 66 + $j ))
phi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com1_${i}_${j},com8_${i}_${j},com13_${i}_${j}


#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 1100 steps per rotation at KAPPA=150 (total torsion 300)
restraint${i}_${j}: ...
    MOVINGRESTRAINT
    ARG=phi${i}_${j},psi${i}_${j} VERSE=B,B
    AT0=-0.035,0.05  STEP0=0        KAPPA0=0,0
    AT1=-0.035,0.05  STEP1=1000     KAPPA1=$KAPPA,$(($KAPPA))
    AT2=-1.5,-1.5  STEP2=1500     KAPPA2=$KAPPA,$(($KAPPA))
    AT3=-3.1,-3.25    STEP3=2100     KAPPA3=$KAPPA,$(($KAPPA))
    AT4=-3.1,-3.25   STEP4=4000     KAPPA4=0,0
...
EOF
done

for ((i=$(( $NMOLS /2)); i<$(($NMOLS/ 2 + $NROT)); i+=1)); do
    factor=$(( $i * $NATOMS + $j ))
        cat >> 'plumed.dat' << EOF
com5_${i}_${j}: COM ATOMS=$((117+$factor)),$((118+$factor))
com1_${i}_${j}: COM ATOMS=$((66+$factor)),$((67+$factor)),$((109+$factor))
com8_${i}_${j}: COM ATOMS=$((63+$factor)),$((64+$factor))
com13_${i}_${j}: COM ATOMS=$((43+$factor)),$((44+$factor)),$((45+$factor)),$((50+$factor))
phi${i}_${j}: TORSION ATOMS=$(( 63 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 110 + $j ))
psi${i}_${j}: TORSION ATOMS=$(( 55 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 66 + $j ))
phi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com1_${i}_${j},com8_${i}_${j},com13_${i}_${j}


#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 1100 steps per rotation at KAPPA=150 (total torsion 300)
restraint${i}_${j}: ...
    MOVINGRESTRAINT
    ARG=phi${i}_${j},psi${i}_${j} VERSE=B,B
    AT0=0.035,-0.05   STEP0=0        KAPPA0=0,0
    AT1=0.035,-0.05 STEP1=1000     KAPPA1=$KAPPA,$(($KAPPA))
    AT2=1.5,1.5   STEP2=1500     KAPPA2=$KAPPA,$(($KAPPA))
    AT3=3.1,3.25    STEP3=2100     KAPPA3=$KAPPA,$(($KAPPA))
    AT4=3.1,3.25    STEP4=4000     KAPPA4=0,0
...
EOF
    done

j=$(( $NMOLS * $NATOMS + 2*$NMOLS ))
   for ((i=0; i<$(( $NROT )); i+=1)); do
    factor=$(( $i * $NATOMS + $j ))
        cat >> 'plumed.dat' << EOF
com5_${i}_${j}: COM ATOMS=$((117+$factor)),$((118+$factor))
com1_${i}_${j}: COM ATOMS=$((66+$factor)),$((67+$factor)),$((109+$factor))
com8_${i}_${j}: COM ATOMS=$((63+$factor)),$((64+$factor))
com13_${i}_${j}: COM ATOMS=$((43+$factor)),$((44+$factor)),$((45+$factor)),$((50+$factor))
phi${i}_${j}: TORSION ATOMS=$(( 63 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 110 + $j ))
psi${i}_${j}: TORSION ATOMS=$(( 55 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 66 + $j ))
phi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com1_${i}_${j},com8_${i}_${j},com13_${i}_${j}


#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 1100 steps per rotation at KAPPA=150 (total torsion 300)
restraint${i}_${j}: ...
    MOVINGRESTRAINT
    ARG=phi${i}_${j},psi${i}_${j} VERSE=B,B
    AT0=-0.035,0.05   STEP0=0        KAPPA0=0,0
    AT1=-0.035,0.05   STEP1=1000     KAPPA1=$KAPPA,$(($KAPPA))
    AT2=-1.5,-1.5   STEP2=1500     KAPPA2=$KAPPA,$(($KAPPA))
    AT3=-3.1,-3.25    STEP3=2100     KAPPA3=$KAPPA,$(($KAPPA))
    AT4=-3.1,-3.25    STEP4=4000    KAPPA4=0,0
...
EOF
done

for ((i=$(( $NMOLS /2)); i<$(($NMOLS/ 2 + $NROT)); i+=1)); do
    factor=$(( $i * $NATOMS + $j ))
        cat >> 'plumed.dat' << EOF
com5_${i}_${j}: COM ATOMS=$((117+$factor)),$((118+$factor))
com1_${i}_${j}: COM ATOMS=$((66+$factor)),$((67+$factor)),$((109+$factor))
com8_${i}_${j}: COM ATOMS=$((63+$factor)),$((64+$factor))
com13_${i}_${j}: COM ATOMS=$((43+$factor)),$((44+$factor)),$((45+$factor)),$((50+$factor))
phi${i}_${j}: TORSION ATOMS=$(( 63 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 110 + $j ))
psi${i}_${j}: TORSION ATOMS=$(( 55 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 66 + $j ))
phi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com1_${i}_${j},com8_${i}_${j},com13_${i}_${j}


#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 1100 steps per rotation at KAPPA=150 (total torsion 300)
restraint${i}_${j}: ...
    MOVINGRESTRAINT
    ARG=phi${i}_${j},psi${i}_${j} VERSE=B,B
    AT0=0.035,-0.05   STEP0=0        KAPPA0=0,0
    AT1=0.035,-0.05   STEP1=1000     KAPPA1=$KAPPA,$(($KAPPA))
    AT2=1.5,1.5   STEP2=1500     KAPPA2=$KAPPA,$(($KAPPA))
    AT3=3.1,3.25    STEP3=2100     KAPPA3=$KAPPA,$(($KAPPA))
    AT4=3.1,3.25    STEP4=4000     KAPPA4=0,0
...
EOF
    done


printf "PRINT STRIDE=10000 FILE=../COLVARS/COLVAR_${KAPPA}_${tropic} ARG=*" >> plumed.dat


#Save rotation args
gmx grompp -f mdrun_fix.mdp -c ../structure_files/fibre_run.gro -t ../structure_files/fibre_run.cpt \
    -p ../run/topol.top  -o ../run/mdrun${KAPPA}_${tropic}.tpr -maxwarn 2
gmx mdrun -deffnm ../run/mdrun${KAPPA}_${tropic} -plumed plumed.dat -v -nt $NCORES -gpu_id $GPU_ID
printf "0\n" | gmx editconf -f ../run/mdrun${KAPPA}_${tropic}.gro -o ../run/mdrun${KAPPA}_${tropic}.pdb -n ../structure_files/fibre.ndx 