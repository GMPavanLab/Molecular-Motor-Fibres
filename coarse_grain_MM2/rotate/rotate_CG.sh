#!/bin/bash
#Not up to date
NCORES=$1
GPU_ID=$2
KAPPA=$3
MOL=MM2

mkdir -p CG
rm -f CG/*
cp ../match-dist/v5/run/run.gro CG/MM2.gro
cp ../match-dist/v5/run/run.cpt CG/MM2.cpt
cp ../match-dist/v5/CG_system.top CG/topol.top
cp ../match-dist/v5/CG_MM2_v5.itp CG/CG_MM2_v5.itp 

cd CG
#Steered MD
cat > plumed.dat << EOF
dist_A: RMSD REFERENCE=../../structure_files/CG_MM2a_lit.pdb TYPE=OPTIMAL
dist_B: RMSD REFERENCE=../../structure_files/CG_MM2b_lit.pdb TYPE=OPTIMAL 
dist_I1: RMSD REFERENCE=../../structure_files/CG_MM2I1_lit.pdb TYPE=OPTIMAL

phi: TORSION ATOMS=5,7,8,11
psi: TORSION ATOMS=5,1,8,13
theta: TORSION ATOMS=5,7,8,9

#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint: MOVINGRESTRAINT ...
    ARG=phi,psi,theta VERSE=B,B,B
    AT0=2.1,-0.7,pi   STEP0=1000        KAPPA0=0,0,0
    AT1=-2.45,-1.4,3   STEP1=1500     KAPPA1=0,$KAPPA,100
    AT2=-2.45,-1.75,pi   STEP2=2000     KAPPA2=$KAPPA,$KAPPA,100
    AT3=-2.45,-1.75,pi   STEP3=5000     KAPPA3=$KAPPA,$KAPPA,100
    AT4=-2.45,-1.75,pi   STEP4=10000     KAPPA4=0,0,0
...
PRINT STRIDE=1 ARG=phi,psi,theta,dist_A,dist_I1,dist_B FILE=COLVAR
EOF

gmx grompp -f ../mdp/run-CG.mdp -c MM2.gro -t MM2.cpt -p topol.top -o run-CG.tpr -maxwarn 1
gmx mdrun -deffnm run-CG -nt $NCORES -gpu_id $GPU_ID -v -plumed plumed.dat
cd ..
