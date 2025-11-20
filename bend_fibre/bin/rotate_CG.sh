#!/bin/bash
#WORKS USING RANDOM SEED  -1512343793
LDSEED=-1
while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--file_name)
      shift #past argument
      FILE_PATH=$1
      shift
      ;;
    -t|--checkpoint)
        shift #past argument
        CHECKPOINT=$1
        shift
        ;;
    -n|--index)
        shift #past argument
        INDEX=$1
        shift
        ;;
    -p|--topology)
      shift #past argument
      TOPOLOGY_PATH=$1
      shift
      ;;
    -deffnm|--define_name)
       NAME=$2
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
    -chunks)
      start=$2
      end=$3
      shift;shift;shift
      ;;
    -seed)
      LDSEED=$2
      shift;shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      shift;
      ;;
    *)
      echo "Unknown value $1"
      return 1
      shift;
      ;;
  esac
done

MOL=MM2
NMOLS=440
NATOMS=28
KAPPA=65
STEPLEN=5000
START=1000

mkdir -p ../run
cp $TOPOLOGY_PATH ../run/${NAME}_rotate_topol.top

rm ../run/run/\#*
rm step*
mkdir -p ../output/COLVARS
rm -f ../output/COLVARS/${NAME}*


#Measure state of one molecule MD
cat > plumed.dat << EOF
#Steered MD
EOF

#Steered MD
for j in 0 55 110 165 220 275 330 385; do
    for ((i=$start; i<$end; i+=1)); do
    factor=$(( $i * $NATOMS + $j*$NATOMS ))
    RSp=""
    RSm="-"
    RSL="L"
    #If j == 55 or 165, RSp = "-" RSm = ""
    if [ $j -eq 55 ] || [ $j -eq 165 ] || [ $j -eq 275 ] || [ $j -eq 385 ]; then
        RSp="-"
        RSm=""
        RSL="U"
    fi

        cat >> 'plumed.dat' << EOF
alpha${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor + 7 )) AXIS=$(( $factor + 1 )),$(( $factor + 3 )) VECTOR2=$(( $factor + 8 )),$(( $factor + 14 ))
beta${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 1 )),$(( $factor + 14 ))
alpha2_${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +7 )) AXIS=$(( $factor + 1 )),$(( $factor +3 )) VECTOR2=$(( $factor + 8 )),$(( $factor + 13 ))
beta2_${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 1 )),$(( $factor + 13 ))
gamma${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 1 )),$(( $factor + 8 ))
gamma2_${i}_${j}: TORSION VECTOR1=$(( $factor + 1 )),$(( $factor +3 )) AXIS=$(( $factor + 1 )),$(( $factor + 7 )) VECTOR2=$(( $factor + 13 )),$(( $factor + 14 ))


#Rotation speed in article: 80 deg/ps -> 1.3962634 rad/ps -> 2.17 / 0.002 = 110 steps per rotation at KAPPA=300
restraint${i}_${j}: MOVINGRESTRAINT ...
    ARG=alpha${i}_${j},alpha2_${i}_${j},beta${i}_${j},beta2_${i}_${j},gamma${i}_${j},gamma2_${i}_${j} VERSE=B,B,B,${RSL},${RSL},B
    AT0=${RSp}0.5,${RSp}0.5,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP0=$START                   KAPPA0=0,0,0,0,0,0
    AT1=${RSp}0.5,${RSp}0.5,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP1=$(($START + 1000))       KAPPA1=$KAPPA,$KAPPA,$KAPPA,$KAPPA,30,50
    AT2=${RSp}2.8,${RSp}2.8,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP2=$(($START + 1000 + 2*$STEPLEN))       KAPPA2=$KAPPA,$KAPPA,$KAPPA,$KAPPA,30,50
    AT3=${RSp}2.8,${RSp}2.8,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP3=$((35*$STEPLEN))       KAPPA3=$KAPPA,$KAPPA,$KAPPA,$KAPPA,30,50
    AT4=${RSp}2.8,${RSp}2.8,${RSm}2,${RSm}2,${RSm}2.8,${RSp}pi   STEP4=$((36*$STEPLEN))       KAPPA4=0,0,0,0,0,0
...
EOF
done
done
#AT1=${RSm}2.45,${RSm}1.4,${RSp}0.9  STEP1=$((3000 + 3*$STEPLEN )) KAPPA1=$KAPPA,$KAPPA,$KAPPA
#AT2=${RSm}2.45,${RSm}1.4,${RSp}0.9  STEP2=$((3000 + 4*$STEPLEN )) KAPPA2=$KAPPA,$KAPPA,$KAPPA
#Compute radius
for ((i=0; i<$NMOLS; i+=1)); do
    cat >> 'plumed.dat' << EOF

hed$i: CENTER ATOMS=$(( $i*$NATOMS + 25 ))-$(( $i*$NATOMS + 28 ))
core$i: CENTER ATOMS=$(( $i*$NATOMS + 1 ))-$(( $i*$NATOMS + 14 ))
tail$i: CENTER ATOMS=$(( $i*$NATOMS + 15 ))-$(( $i*$NATOMS + 24 ))
oxy$i: CENTER ATOMS=$(( 17 + $i * $NATOMS )),$(( 18 + $i * $NATOMS ))
EOF
done

#Dump positions for radii calculation
for name in hed core tail oxy; do
    printf "DUMPATOMS STRIDE=100 FILE=../output/COLVARS/${NAME}_${name}.xyz ATOMS=${name}0" >> plumed.dat
    for ((i=1; i<$(($NMOLS)); i+=1)); do
        printf ",${name}$i" >> plumed.dat
    done
    printf "\n" >> plumed.dat
done

NEDIT=8
diff=$(($end - $start))
if ! grep -q "MOTORB-R" ../run/${NAME}_rotate_topol.top; then
  count=0

  while IFS= read -r line; do
    if [[ "$line" =~ MOTOR-R ]]; then
      count=$((count+1))
      if [ $count -le ${NEDIT} ]; then
        echo "MOTORB-R ${diff}" >> ../tmp/tmp.top
      fi
    fi
    if [[ "$line" =~ MOTOR-S ]]; then
      count=$((count+1))
      if [ $count -le ${NEDIT} ]; then
        echo "MOTORB-S ${diff}" >> ../tmp/tmp.top
      fi
    fi
    echo $line >> ../tmp/tmp.top
  done <  ../run/${NAME}_rotate_topol.top

  mv ../tmp/tmp.top  ../run/${NAME}_rotate_topol.top
  awk 'BEGIN {count=0} /MOTOR-[RS]/ { if (count < '"${NEDIT}"') { $2 = $2 - '"$diff"'; count++ } } 1' ../run/${NAME}_rotate_topol.top > ../tmp/tmp.top
else
awk 'BEGIN {count=0} /MOTOR-[RS]/ { if (count < '"${NEDIT}"') { $2 = $2 + '"$diff"'; count++ } } 1' ../run/${NAME}_rotate_topol.top > ../tmp/tmp2.top
awk 'BEGIN {count=0} /MOTOR-[RS]/ { if (count < '"${NEDIT}"') { $2 = $2 - '"$diff"'; count++ } } 1' ../tmp/tmp2.top > ../tmp/tmp.top
fi
mv ../tmp/tmp.top ../run/${NAME}_rotate_topol.top

#Run MD
mkdir -p ../run/run/
sed -i -E '/ld-seed/s/(-?[0-9]+)/'"$LDSEED"'/' mdp/run.mdp
gmx grompp -p ../run/${NAME}_rotate_topol.top -n $INDEX -f mdp/run.mdp -c $FILE_PATH -t $CHECKPOINT -o ../run/run/${NAME}.tpr -maxwarn 1
gmx mdrun -nt $NCORES -gpu_id $GPU_ID -deffnm ../run/run/${NAME} -plumed plumed.dat -v
