#!/bin/bash/
: '
Runs a 500 ns metadynamics simulation for a pair of motors, either AA or CG. The metadynamics are placed on the distance between COM to compute dimerisation PMF.
Weights are then reweighted based on the energy profile on the last frame.
The different motors are namned as follows:
    A/B: MM2a/MM2b
    R/S: Isomer
    C: CG

INPUT:
    $1: Number of cores to run on 
    $2: GPU ID to run on
    $3: The pair to run the simulation on. (AR-AR, AR-AS, AR-BR, AR-BS, BR-BR, BR-BS)
OUTPUT:
    -metad/COLVAR_${ON}/: COLVAR files for each simulation
    -metad/(em|eq*|npt|nvt|run): Gromacs output files for each simulation
'
NCORES=$1
GPU_ID=$2
ON=$3

AR=../structure_files/MM2a-R.gro
AS=../structure_files/MM2a-S.gro
BR=../structure_files/MM2b-R.gro
BS=../structure_files/MM2b-S.gro
CGAR=../structure_files/CG_MM2a-R.gro
CGAS=../structure_files/CG_MM2a-S.gro
CGBR=../structure_files/CG_MM2b-R.gro
CGBS=../structure_files/CG_MM2b-S.gro
TOP=../structure_files/MM2.top
CGTOPA=../structure_files/CG_MM2a-R.top
CGTOPB=../structure_files/CG_MM2b-R.top

case $ON in
    AR-AR)
        name=AR-AR
        echo "Running $name"
        bash run_AAmetad.sh $NCORES $GPU_ID $AR $AR $TOP MM2-R $name
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    AR-AS)
        name=AR-AS
        echo "Running $name"
        bash run_AAmetad.sh $NCORES $GPU_ID $AR $AS $TOP MM2-S $name
        mkdir -p ../metad/COLVARS_${name}
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    AR-BR)
        name=AR-BR
        echo "Running $name"
        bash run_AAmetad.sh $NCORES $GPU_ID $AR $BR $TOP MM2-R $name
        mkdir -p ../metad/COLVARS_${name}
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    AR-BS)
        name=AR-BS
        echo "Running $name"
        bash run_AAmetad.sh $NCORES $GPU_ID $AR $BS $TOP MM2-S $name
        mkdir -p ../metad/COLVARS_${name}
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    BR-BR)
        name=BR-BR
        echo "Running $name"
        bash run_AAmetad.sh $NCORES $GPU_ID $BR $BR $TOP MM2-R $name
        mkdir -p ../metad/COLVARS_${name}
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    BR-BS)
        name=BR-BS
        echo "Running $name"
        bash run_AAmetad.sh $NCORES $GPU_ID $BR $BS $TOP MM2-S $name
        mkdir -p ../metad/COLVARS_${name}
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    CG_AR-AR)
        name=CG_AR-AR
        echo "Running $name"
        bash run_CGmetad.sh $NCORES $GPU_ID $CGAR $CGAR $CGTOPA MM2a-R $name
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    CG_AR-AS)
        name=CG_AR-AS
        echo "Running $name"
        bash run_CGmetad.sh $NCORES $GPU_ID $CGAR $CGAS $CGTOPA MM2a-S $name
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    CG_AR-BR)
        name=CG_AR-BR
        echo "Running $name"
        bash run_CGmetad.sh $NCORES $GPU_ID $CGAR $CGBR $CGTOPA MM2b-R $name
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    CG_BR-BR)
        name=CG_BR-BR
        echo "Running $name"
        bash run_CGmetad.sh $NCORES $GPU_ID $CGBR $CGBR $CGTOPB MM2b-R $name
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    CG_BR-BS)
        name=CG_BR-BS
        echo "Running $name"
        bash run_CGmetad.sh $NCORES $GPU_ID $CGBR $CGBS $CGTOPB MM2b-S $name
        bash compute_potentials.sh ../metad/run/md_${name}.xtc ../metad/COLVARS_${name}/HILLS $name
        echo "Done"
        ;;
    esac
rm ../metad/*/\#*
rm ../metad/*/bck*