#!/bin/bash
#Not up to date
NCORES=$1
GPU_ID=$2
MOL=MM2

mkdir -p AA
rm AA/*
cp ../AA-reference/npt/npt_MM2.gro AA/MM2.gro
cp ../AA-reference/npt/npt_MM2.cpt AA/MM2.cpt
cp AA.ndx AA/index.ndx
cp MM2-AA.top AA/topol.top


cd AA
#Steered MD
cat > plumed.dat << EOF
phi: TORSION ATOMS=25,26,27,39
psi: TORSION ATOMS=17,26,27,28

#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 1100 steps per rotation at KAPPA=300
restraint: ...
    MOVINGRESTRAINT
    ARG=phi
    AT0=-0.035    STEP0=0        KAPPA0=0
    AT1=-0.035    STEP1=1000     KAPPA1=300
    AT2=-3.063    STEP2=2100     KAPPA2=300
    AT3=-3.063    STEP3=3000     KAPPA3=0
...
PRINT STRIDE=10 ARG=* FILE=COLVAR
EOF

gmx grompp -f ../mdp/run-AA.mdp -c MM2.gro -t MM2.cpt -p topol.top -o run-AA.tpr -n index.ndx -maxwarn 1
gmx mdrun -deffnm run-AA -nt $NCORES -gpu_id $GPU_ID -v -plumed plumed.dat
