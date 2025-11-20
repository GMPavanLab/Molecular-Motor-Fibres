#!/bin/bash

mkdir -p ../run
cp ../structure_files/fibre.top ../run/topol.top
rm COLVAR core* hed* mdout.mdp oxy* plumed* bck* step*
rm ../run/\#*
NCORES=$1
GPU_ID=$2
KAPPA=$3
tropic=single
NMOLS=110
NATOMS=152
cat > 'plumed.dat' << EOF
MOLINFO STRUCTURE=../structure_files/fibre_structure.pdb
EOF

#make motors whole
printf 'WHOLEMOLECULES ENTITY0=1-152 ENTITY1=8360-8512\n' >> plumed.dat

#Perform steered metad on dihedrals
for j in 0; do
    for i in 0; do
    factor=$(( $i * $NATOMS + $j ))
        cat >> 'plumed.dat' << EOF
com5_${i}_${j}: COM ATOMS=$((117+$factor)),$((118+$factor))
com1_${i}_${j}: COM ATOMS=$((66+$factor)),$((67+$factor)),$((109+$factor))
com7_${i}_${j}: COM ATOMS=$((110+$factor)),$((111+$factor)),$((116+$factor))
com8_${i}_${j}: COM ATOMS=$((63+$factor)),$((64+$factor))
com11_${i}_${j}: COM ATOMS=$((39+$factor)),$((40+$factor)),$((51+$factor))
com13_${i}_${j}: COM ATOMS=$((43+$factor)),$((44+$factor)),$((45+$factor)),$((50+$factor))

phi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com1_${i}_${j},com8_${i}_${j},com13_${i}_${j}
psi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com7_${i}_${j},com8_${i}_${j},com11_${i}_${j}
phi${i}_${j}: TORSION ATOMS=$(( 63 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 110 + $j ))
psi${i}_${j}: TORSION ATOMS=$(( 55 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 66 + $j ))

#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint${i}_${j}: ...
    MOVINGRESTRAINT
    ARG=phi${i}_${j},psi${i}_${j},phi_CG${i}_${j} VERSE=B,B,L
    AT0=-0.035,0.05,-1.2    STEP0=1000        KAPPA0=0,0,600
    AT1=-0.035,0.05,-1.2   STEP1=2000     KAPPA1=$KAPPA,$KAPPA,600
    AT2=-1.5,-1.5,-1.2    STEP2=2500     KAPPA2=$KAPPA,$KAPPA,600
    AT3=-3.1,-3.25,-1.2    STEP3=3000     KAPPA3=$KAPPA,$KAPPA,600
    AT4=-3.1,-3.25,-1.2    STEP4=4000     KAPPA4=0,0,0
...
EOF
done

for i in 55; do
    factor=$(( $i * $NATOMS + $j ))
cat >> 'plumed.dat' << EOF

com5_${i}_${j}: COM ATOMS=$((117+$factor)),$((118+$factor))
com1_${i}_${j}: COM ATOMS=$((66+$factor)),$((67+$factor)),$((109+$factor))
com7_${i}_${j}: COM ATOMS=$((110+$factor)),$((111+$factor)),$((116+$factor))
com8_${i}_${j}: COM ATOMS=$((63+$factor)),$((64+$factor))
com11_${i}_${j}: COM ATOMS=$((39+$factor)),$((40+$factor)),$((51+$factor))
com13_${i}_${j}: COM ATOMS=$((43+$factor)),$((44+$factor)),$((45+$factor)),$((50+$factor))

phi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com1_${i}_${j},com8_${i}_${j},com13_${i}_${j}
psi_CG${i}_${j}: TORSION ATOMS=com5_${i}_${j},com7_${i}_${j},com8_${i}_${j},com11_${i}_${j}

phi${i}_${j}: TORSION ATOMS=$(( 63 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 110 + $j ))
psi${i}_${j}: TORSION ATOMS=$(( 55 + $i * $NATOMS + $j )),$(($i * $NATOMS + 64 + $j )),$(($i * $NATOMS + 65 + $j )),$(($i * $NATOMS + 66 + $j ))

#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint${i}_${j}: ...
    MOVINGRESTRAINT
    ARG=phi${i}_${j},psi${i}_${j},phi_CG${i}_${j} VERSE=B,B,U
    AT0=0.035,-0.05,1.2    STEP0=1000        KAPPA0=0,0,600
    AT1=0.035,-0.05,1.2   STEP1=2000     KAPPA1=$KAPPA,$KAPPA,600
    AT2=1.5,1.5,1.2    STEP2=2500     KAPPA2=$KAPPA,$KAPPA,600
    AT3=3.1,3.25,1.2    STEP3=3000     KAPPA3=$KAPPA,$KAPPA,600
    AT4=3.1,3.25,1.2    STEP4=4000     KAPPA4=0,0,0
...
EOF
    done
done
#Save rotation args
printf "PRINT STRIDE=10 FILE=COLVAR_${KAPPA}_${tropic} ARG=*" >> plumed.dat



gmx grompp -f mdrun.mdp -c ../structure_files/fibre_run.gro -t ../structure_files/fibre_run.cpt \
    -p ../run/topol.top -n ../structure_files/temp_groups.ndx -o ../run/mdrun${KAPPA}_${tropic}.tpr -maxwarn 2
gmx mdrun -deffnm ../run/mdrun${KAPPA}_${tropic} -plumed plumed.dat -v -nt $NCORES -nb cpu -gpu_id $GPU_ID