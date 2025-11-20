#!/bin/bash
: '
Reweights a metadynamics simulation using the bias of the last frame. 
INPUT:
    $1: .xtc - file
        Path to the trajectory file of the metadynamics simulation
    $2: path
        Path to the HILLS file of the metadynamics simulation
    $3: str
        Name of the simulation. Must be the same as the name used when running metadynamics.
'
traj=$1
HILLS=$2
name=$3
N=49 # Number of atoms in each motor

cp $traj ../metad/COLVARS_${name}/traj.xtc
cp $HILLS ../metad/COLVARS_${name}/HILLS_TO_BE_REWEIGHTED
cd ../metad/COLVARS_${name}

rm COLVAR-${name} dist-${name}.weight
cat << EOF > plumed.dat
RESTART
WHOLEMOLECULES ENTITY0=1-$N, ENTITY1=$(($N+1))-$(($N*2))

com1: CENTER ATOMS=1-$N
com2: CENTER ATOMS=$(($N+1))-$(($N*2))

dist: DISTANCE ATOMS=com1,com2

metad: METAD ARG=dist PACE=100000000 HEIGHT=0.2 SIGMA=0.05 ...
    FILE=HILLS_TO_BE_REWEIGHTED BIASFACTOR=10 TEMP=300 
    GRID_MIN=0.0 GRID_MAX=4.0 GRID_BIN=400
    ...

PRINT STRIDE=1 ARG=dist,metad.bias FILE=COLVAR-${name}
EOF
plumed driver --mf_xtc traj.xtc --plumed plumed.dat

bmax=`awk 'BEGIN{max=0.}{if($1!="#!" && $3>max)max=$3}END{print max}' COLVAR-${name}`

# print phi values and weights
awk '{if($1!="#!") print $2,exp(($3-bmax)/kbt)}' kbt=2.494339 bmax=$bmax COLVAR-${name} > dist-${name}.weight
cd ../../bin