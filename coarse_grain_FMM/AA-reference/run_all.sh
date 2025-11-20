#!/bin/bash
: '
Creates AA-references for coarse graining by running a 1000 ns simulation.
INPUT:
    $1: int
        Number of cores to use
    $2: int
        GPU ID to use
    $3: str
        Name of the molecule to simulate. ../structure_files/${MOL}.(gro|top) must both exist
OUTPUT:
    -run/* : Gromacs output files for the simulation
'

NCORES=16
GPU_ID=0
iso=R
NSTEPS=-2
while [[ $# -gt 0 ]]; do
  case $1 in
    -nc|--n_cores) #Path to mol2 file
      shift #past argument
      NCORES=$1
      shift
      ;;
    -gpu_id) #What the files should be named
      shift #past argument
      GPU_ID=$1
      shift
      ;;
    -iso) #What the files should be named
      shift #past argument
      iso=$1
      shift
      ;;
    -ns|--n_steps) #Path to main_chain file
      shift #past argument
      NSTEPS=$1
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
  esac
done
mkdir -p tmp/
mkdir -p logs/
mkdir -p output/

cp ../structure_files/motor2m-${iso}.top topol2m-${iso}.top
echo "RUNNING MOTOR2m-${iso}"
echo "Editing conf"
gmx editconf -f ../structure_files/motor2m-${iso}.gro -o tmp/motor2m-${iso}_box.gro -c -d 1.0 -bt dodecahedron > logs/edit_conf2m-${iso}.out

TOPOL=topol2m-${iso}.top

gmx solvate -cp tmp/motor2m-${iso}_box.gro -cs spc216.gro -o tmp/motor2m-${iso}_solv.gro -p $TOPOL > logs/solvate2m-${iso}.out

echo "Adding ions"
gmx grompp -f mdp/ions.mdp -c tmp/motor2m-${iso}_solv.gro -p $TOPOL -o tmp/ions.tpr
printf '6\n0\n' | gmx genion -s tmp/ions.tpr -o tmp/motor2m-${iso}_ions.gro -p $TOPOL -pname NA -nname CL -neutral

echo "Making index"
printf 'name 1 MOTOR\n7 | 9\nq\n' | gmx make_ndx -f tmp/motor2m-${iso}_ions.gro -o ../structure_files/index-${iso}.ndx
INDEX=../structure_files/index-${iso}.ndx

mkdir -p em
echo "Energy minimisation"
gmx grompp -f mdp/emin.mdp -c tmp/motor2m-${iso}_ions.gro -p $TOPOL -o em/em2m-${iso}.tpr > logs/em2m-${iso}_grompp.out
gmx mdrun -deffnm em/em2m-${iso}  -nt $NCORES > logs/em2m-${iso}_run.out 

printf '10\n0\n' | gmx energy -f em/em2m-${iso}.edr -o output/em2m-${iso}_potential.xvg > logs/em2m-${iso}_energy.out
echo "Done, potential in logs"

mkdir -p nvt
echo "NVT"
gmx grompp -f mdp/nvt.mdp -c em/em2m-${iso}.gro -r em/em2m-${iso}.gro -p $TOPOL -n $INDEX -o nvt/nvt2m-${iso}.tpr > logs/nvt2m-${iso}_grompp.out
gmx mdrun -deffnm nvt/nvt2m-${iso} -nt $NCORES -gpu_id $GPU_ID > logs/nvt2m-${iso}_run.out

printf '15\n0\n' | gmx energy -f nvt/nvt2m-${iso}.edr -o output/nvt2m-${iso}_temperature.xvg > logs/nvt2m-${iso}_energy.out
echo "Done, temperature in logs"

mkdir -p npt
echo "NPT"
gmx grompp -f mdp/npt.mdp -c nvt/nvt2m-${iso}.gro -r nvt/nvt2m-${iso}.gro -t nvt/nvt2m-${iso}.cpt -p $TOPOL -n $INDEX -o npt/npt2m-${iso}.tpr > logs/npt2m-${iso}_grompp.out
gmx mdrun -deffnm npt/npt2m-${iso} -nt $NCORES -gpu_id $GPU_ID > logs/npt2m-${iso}_run.out

printf '17\n0\n' | gmx energy -f npt/npt2m-${iso}.edr -o output/npt2m-${iso}_pressure.xvg > logs/npt2m-${iso}_energy.out
printf '23\n0\n' | gmx energy -f npt/npt2m-${iso}.edr -o output/npt2m-${iso}_density.xvg > logs/npt2m-${iso}_energy2.out

echo "Done, pressure and density in logs"

mkdir -p run
echo "MD Run 1 ns"
gmx grompp -f mdp/run.mdp -c npt/npt2m-${iso}.gro -t npt/npt2m-${iso}.cpt -p $TOPOL -n $INDEX -o run/md2m-${iso}_run.tpr > logs/md2m-${iso}_grompp.out

gmx mdrun -s run/md2m-${iso}_run.tpr -o run/md2m-${iso}.cpt -c run/md2m-${iso}.gro -x run/md2m-${iso}.xtc -e run/md2m-${iso}.edr \
    -nt $NCORES -gpu_id $GPU_ID -nsteps $NSTEPS> logs/md2m-${iso}_run.out
echo "Done"
rm -r \#*
