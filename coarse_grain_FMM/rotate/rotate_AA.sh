#!/bin/bash

NCORES=$1
GPU_ID=$2
MOL=MM2
KAPPA=300

mkdir -p AA
rm -f AA/*
cp ../AA-reference/npt/npt2m-R.gro AA/FMM.gro
cp ../AA-reference/npt/npt2m-R.cpt AA/FMM.cpt
cp ../structure_files/index-R.ndx AA/index.ndx


cd AA
#Steered MD
cat > plumed.dat << EOF
com1: COM ATOMS=66,67,109
com3: COM ATOMS=70,71
com4: COM ATOMS=65,72
com5: COM ATOMS=117,118
com7: COM ATOMS=110,111,116
com8: COM ATOMS=63,64
com9: COM ATOMS=52,53
com10: COM ATOMS=54,55,56,57,58,59,60,61,62
com11: COM ATOMS=39,40,51
com13: COM ATOMS=43,44,45,50

phi_CG: TORSION ATOMS=com5,com1,com8,com13
psi_CG: TORSION ATOMS=com5,com7,com8,com11
theta: TORSION ATOMS=56,64,65,111
phi: TORSION ATOMS=63,64,65,110
psi: TORSION ATOMS=55,64,65,66
theta_CG: TORSION ATOMS=com5,com1,com8,com10

#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint: MOVINGRESTRAINT ...
    ARG=psi_CG,theta_CG,phi_CG VERSE=B,U,L
    AT0=2,-3,-1.2   STEP0=0        KAPPA0=0,0,600
    AT1=2,2,-1.2   STEP1=1000     KAPPA1=$KAPPA,$KAPPA,600
    AT2=0.5,0.8,-1.2   STEP2=1500     KAPPA2=$KAPPA,$KAPPA,600
    AT3=-2.5,0.8,-1.2    STEP3=2000     KAPPA3=$KAPPA,$KAPPA,600
    AT4=-2.5,0.8,-1.2    STEP4=3000     KAPPA4=0,0,0
...

PRINT STRIDE=1 ARG=phi,psi,phi_CG,psi_CG,theta,theta_CG,alpha_CG FILE=COLVAR
EOF

gmx grompp -f ../mdp/run-AA.mdp -c FMM.gro -t FMM.cpt -p ../topol2m-R.top -o run-AA.tpr -n index.ndx -maxwarn 1
gmx mdrun -deffnm run-AA -nt $NCORES -gpu_id $GPU_ID -v -plumed plumed.dat
echo 1 | gmx editconf -f run-AA.gro -o run-AA.pdb -n index.ndx
cd ..

:'
restraint: MOVINGRESTRAINT ...
    ARG=phi,psi,phi_CG VERSE=B,B,L
    AT0=-0.035,0.05,-1.2   STEP0=0        KAPPA0=0,0,600
    AT1=-0.035,0.05,-1.2   STEP1=1000     KAPPA1=$KAPPA,$KAPPA,600
    AT2=-1.5,-1.5,-1.2   STEP2=1500     KAPPA2=$KAPPA,$KAPPA,600
    AT3=-3.1,-3.25,-1.2    STEP3=2000     KAPPA3=$KAPPA,$KAPPA,600
    AT4=-3.1,-3.25,-1.2    STEP4=3000     KAPPA4=0,0,0
...'