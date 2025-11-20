#!/bin/bash/

ENERGY=$1
OUTPATH=$2

echo 15 | gmx energy -f $ENERGY -o ${OUTPATH}_x.xvg
echo 16 | gmx energy -f $ENERGY -o ${OUTPATH}_y.xvg
echo 17 | gmx energy -f $ENERGY -o ${OUTPATH}_z.xvg
echo 18 | gmx energy -f $ENERGY -o ${OUTPATH}_volume.xvg
echo 19 | gmx energy -f $ENERGY -o ${OUTPATH}_density.xvg