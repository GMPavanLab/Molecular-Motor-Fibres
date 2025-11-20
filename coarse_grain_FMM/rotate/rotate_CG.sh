#!/bin/bash
#-16793629

NCORES=$1
GPU_ID=$2
KAPPA=65
STEPLEN=5000
MOL=MM2
START=5000
mkdir -p CG
rm -f CG/*
cp ../match-dist/v4/run/run.gro CG/FMM.gro
cp ../match-dist/v4/run/run.cpt CG/FMM.cpt
cp ../match-dist/v4/index.ndx CG/index.ndx
cp ../match-dist/v4/CG_system.top CG/topol.top

cp ../match-dist/v4/CG_motor2m-R_v4.itp CG/CG_motor2m-R_v4.itp 

cd CG
#Steered MD
cat > plumed.dat << EOF
dist_A: RMSD REFERENCE=../../structure_files/CG_MM2a_lit.pdb TYPE=OPTIMAL
dist_B: RMSD REFERENCE=../../structure_files/CG_MM2b_lit.pdb TYPE=OPTIMAL 
dist_I1: RMSD REFERENCE=../../structure_files/CG_MM2I1_lit.pdb TYPE=OPTIMAL


#phi: TORSION ATOMS=5,7,8,11
#psi: TORSION ATOMS=5,1,8,13
#omega: TORSION VECTOR1=1,7 AXIS=1,3 VECTOR2=8,10
#theta: TORSION ATOMS=5,7,8,9
alpha: TORSION VECTOR1=1,7 AXIS=1,3 VECTOR2=8,14
beta: TORSION VECTOR1=1,3 AXIS=1,7 VECTOR2=1,14
alpha2: TORSION VECTOR1=1,7 AXIS=1,3 VECTOR2=8,13
beta2: TORSION VECTOR1=1,3 AXIS=1,7 VECTOR2=1,13
#beta3: TORSION VECTOR1=1,3 AXIS=1,7 VECTOR2=11,14
gamma: TORSION VECTOR1=1,3 AXIS=1,7 VECTOR2=1,8
gamma2: TORSION VECTOR1=1,3 AXIS=1,7 VECTOR2=13,14

#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint: MOVINGRESTRAINT ...
    ARG=alpha,alpha2,beta,beta2,gamma,gamma2 VERSE=B,B,B,L,L,B
    AT0=0.5,0.5,-2,-2,-2.8,pi  STEP0=$START        KAPPA0=0,0,0,0,0,0
    AT1=0.5,0.5,-2,-2,-2.8,pi   STEP1=$(($START + 1000))     KAPPA1=$KAPPA,$KAPPA,$KAPPA,$KAPPA,30,50
    AT2=2.8,2.8,-2,-2,-2.8,pi   STEP2=$(($START + 1000 + 2*$STEPLEN )) KAPPA2=$KAPPA,$KAPPA,$KAPPA,$KAPPA,30,50
    AT3=2.8,2.8,-2,-2,-2.8,pi   STEP3=$(($START + 1000 + 3*$STEPLEN )) KAPPA3=$KAPPA,$KAPPA,$KAPPA,$KAPPA,30,50
...


PRINT STRIDE=1 ARG=gamma,gamma2,alpha,beta,alpha2,beta2,dist_A,dist_I1,dist_B FILE=COLVAR
EOF

sed -i '/; Include chain topologies/a #include "../../structure_files/CG_motorb-R_v4.itp"' topol.top
sed -i 's/MOTOR-R/MOTORB-R/g' topol.top

gmx grompp -f ../mdp/run-CG.mdp -c FMM.gro -t FMM.cpt -p topol.top -o run-CG.tpr -n index.ndx -maxwarn 1
gmx mdrun -deffnm run-CG -nt $NCORES -gpu_id $GPU_ID -plumed plumed.dat 
echo 1 | gmx editconf -f run-CG.gro -o run-CG.pdb -n index.ndx
cd ..

:'
restraint: MOVINGRESTRAINT ...
    ARG=phi,psi,theta VERSE=B,B,B
    AT0=2.1,-0.7,pi   STEP0=1000        KAPPA0=0,0,0
    AT1=0,-1.4,3   STEP1=1500     KAPPA1=0,$KAPPA,100
    AT2=-2.45,-1.75,pi   STEP2=2000     KAPPA2=$KAPPA,$KAPPA,100
    AT3=-2.45,-1.75,pi   STEP3=10000     KAPPA3=$KAPPA,$KAPPA,100
    AT4=-2.45,-1.75,pi   STEP4=11000     KAPPA4=0,0,0
...
'