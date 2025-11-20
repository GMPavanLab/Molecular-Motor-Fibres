#!/bin/bash

traj=$1
HILLS=$2
name=$3
N=$4

mkdir -p block_analysis
cp $traj block_analysis/traj.xtc
cp $HILLS block_analysis/HILLS
cd block_analysis

rm COLVAR-${name} dist-${name}.weight
cat << EOF > plumed.dat
RESTART
WHOLEMOLECULES ENTITY0=1-$N, ENTITY1=$(($N+1))-$(($N*2))

com1: CENTER ATOMS=1-$N
com2: CENTER ATOMS=$(($N+1))-$(($N*2))

dist: DISTANCE ATOMS=com1,com2

metad: METAD ARG=dist PACE=100000000 HEIGHT=0.2 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 TEMP=300 
    GRID_MIN=0.0 GRID_MAX=4.0 GRID_BIN=400
    ...

PRINT STRIDE=1 ARG=dist,metad.bias FILE=COLVAR-${name}
EOF
plumed driver --mf_xtc traj.xtc --plumed plumed.dat

bmax=`awk 'BEGIN{max=0.}{if($1!="#!" && $3>max)max=$3}END{print max}' COLVAR-${name}`

# print phi values and weights
awk '{if($1!="#!") print $2,exp(($3-bmax)/kbt)}' kbt=2.494339 bmax=$bmax COLVAR-${name} > dist-${name}.weight
