#!/bin/bash

NCORES=16
GPU_ID=0
FILE_PATH=../run/compress2/compress
OUTPUT_NAME=tmp
solvent_box=../structure_files/box_CG_W_eq.gro


while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--file_name)
      shift #past argument
      FILE_PATH=$1
      shift
      ;;
    -p|--topology)
      shift #past argument
      TOPOLOGY_PATH=$1
      shift
      ;;
    -deffnm|--define_name)
       OUTPUT_NAME=$2
       shift;shift
       ;;
    -nt|--number_cores)
      NCORES=$2
      shift;shift
      ;;
    -gpu_id)
      GPU_ID=$2
      shift;shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
    *)
      echo "Unknown value $1"
      return 1
      ;;
  esac
done

rm */\#*
rm \#*
rm ../run/\#*
rm ../run/*/\#*
rm ../tmp/\#*
rm ../system/\#*

mkdir -p ../system/
mkdir -p ../system/${OUTPUT_NAME}/
mkdir -p ../tmp
mkdir -p ../run
cp $TOPOLOGY_PATH ../tmp/${OUTPUT_NAME}_topol.top
TOPOLOGY_PATH=../tmp/${OUTPUT_NAME}_topol.top

#Equilibriate system and couple isotropically
INDEX=../run/${OUTPUT_NAME}.ndx
gmx make_ndx -f $FILE_PATH -o ${INDEX} <<EOF
2 | 3 | 4 | 5
name 6 SOL_ION
name 7 MOTOR
q
EOF
#Extract box lengths
xybox=$(tail -n 1 $FILE_PATH | awk '{print $1}')
zlen=$(tail -n 1 $FILE_PATH | awk '{print $3}')

mkdir -p ../run/em_vacuum
gmx grompp -p ${TOPOLOGY_PATH} -c $FILE_PATH -f mdp/em.mdp  -o ../run/em_vacuum/em.tpr  -maxwarn 1
gmx mdrun -v -deffnm ../run/em_vacuum/em -gpu_id $GPU_ID -nt $NCORES

mkdir -p ../run/eq1_vacuum
gmx grompp -p ${TOPOLOGY_PATH} -c ../run/em_vacuum/em.gro -n $INDEX -f mdp/eq1.mdp  -o ../run/eq1_vacuum/eq.tpr -maxwarn 1
gmx mdrun -deffnm ../run/eq1_vacuum/eq -gpu_id $GPU_ID -nt $NCORES

for i in {2..5}; do
    mkdir -p ../run/eq${i}_vacuum
    gmx grompp -p ${TOPOLOGY_PATH} -c ../run/eq$(( $i - 1 ))_vacuum/eq.gro -f mdp/eq${i}.mdp \
     -n $INDEX -o ../run/eq${i}_vacuum/eq.tpr  -maxwarn 1
    gmx mdrun -v -deffnm ../run/eq${i}_vacuum/eq -gpu_id $GPU_ID -nt $NCORES -nsteps 100000
done

#Make molecules whole and free fibre from pbc 
gmx grompp -f ../bin/mdp/em.mdp -c ../run/eq5_vacuum/eq.gro -p $TOPOLOGY_PATH -o ../tmp/tmp.tpr -maxwarn 1
echo 0 | gmx trjconv -f ../run/eq5_vacuum/eq.gro -s ../tmp/tmp.tpr -o ../tmp/fibre_whole.gro -pbc whole
gmx editconf -f ../tmp/fibre_whole.gro -c -o ../tmp/box.gro -box $xybox $xybox $(echo "$zlen+10" | bc -l)

#Solvate system
gmx solvate -cp ../tmp/box.gro  -cs $solvent_box \
   -o ../tmp/box_sol.gro -p ${TOPOLOGY_PATH}
python remove_water_inside_fibre.py ../tmp/box_sol.gro ../tmp/box_clean.gro $TOPOLOGY_PATH


#Make index and ease position restraints
printf "[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n" > ../structure_files/posres_fibre.itp
for i in {1..28}; do
  echo "$i 1 200 200 200" >> ../structure_files/posres_fibre.itp
done

INDEX=../run/${OUTPUT_NAME}.ndx
gmx make_ndx -f ../tmp/box_clean.gro -o ${INDEX} <<EOF
2 | 3 | 4 | 5
6 | 7
name 8 MOTOR
name 9 SOL_ION
q
EOF

#Equilibriate water
gmx grompp -p ${TOPOLOGY_PATH} -c ../tmp/box_clean.gro -f mdp/em_solvant.mdp  \
 -r ../tmp/box_clean.gro -o ../run/em_solvant/em_W.tpr  -maxwarn 1
gmx mdrun -deffnm ../run/em_solvant/em_W -gpu_id $GPU_ID -nt $NCORES

gmx grompp -p ${TOPOLOGY_PATH} -c ../run/em_solvant/em_W.gro -f mdp/eq_solvant2.mdp  \
   -r ../tmp/box_clean.gro -o ../run/eq_solvant/eq_W.tpr -maxwarn 1 -n $INDEX
gmx mdrun -deffnm ../run/eq_solvant/eq_W -gpu_id $GPU_ID -nt $NCORES

#Equilibrate system
mkdir -p ../run/em
gmx grompp -p ${TOPOLOGY_PATH} -c ../run/eq_solvant/eq_W.gro -f mdp/em.mdp  -o ../run/em/em.tpr  -maxwarn 1
gmx mdrun -v -deffnm ../run/em/em -gpu_id $GPU_ID -nt $NCORES

mkdir -p ../run/eq1
gmx grompp -p ${TOPOLOGY_PATH} -c ../run/em/em.gro -n $INDEX -f mdp/eq1.mdp  -o ../run/eq1/eq.tpr -maxwarn 1
gmx mdrun -deffnm ../run/eq1/eq -gpu_id $GPU_ID -nt $NCORES

for i in {2..5}; do
    mkdir -p ../run/eq${i}
    gmx grompp -p ${TOPOLOGY_PATH} -c ../run/eq$(( $i - 1 ))/eq.gro -f mdp/eq${i}.mdp \
     -n $INDEX -o ../run/eq${i}/eq.tpr  -maxwarn 1
    gmx mdrun -v -deffnm ../run/eq${i}/eq -gpu_id $GPU_ID -nt $NCORES
done

cp ../run/eq5/eq.gro ../system/$OUTPUT_NAME/equilibrated.gro
cp ../run/eq5/eq.cpt ../system/$OUTPUT_NAME/equilibrated.cpt
cp ${TOPOLOGY_PATH} ../system/${OUTPUT_NAME}_topol_sol.top