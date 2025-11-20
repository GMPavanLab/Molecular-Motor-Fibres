#!/bin/bash/

GPU_ID=$1
NT=$2
NAME=$3
cp ../structure_files/CG_system.top ../tmp/system.top

bash packmol_structure_shifted.sh ../structure_files/fibre_structure.pdb ../tmp/system.gro \
     ../tmp/system.top $GPU_ID

bash build_system.sh -f ../tmp/system.gro -p ../tmp/system.top -deffnm $NAME -gpu_id $GPU_ID -nt $NT

bash equilibriate_shifted.sh -f ../system/${NAME}/compressed.gro -p ../system/${NAME}_topol.top -deffnm $NAME -gpu_id $GPU_ID -nt $NT

bash rotate_CG_shifted.sh -f ../system/$NAME/equilibrated.gro -p ../system/${NAME}_topol_sol.top \
    -t ../system/$NAME/equilibrated.cpt --index ../run/${NAME}.ndx -deffnm ${NAME}5 -gpu_id $GPU_ID -nt $NT \
    -chunks 0 5

bash rotate_CG_shifted.sh -f ../run/run/${NAME}5.gro -p ../run/${NAME}5_rotate_topol.top \
     -t ../run/run/${NAME}5.cpt --index ../run/${NAME}.ndx -deffnm ${NAME}25 -gpu_id $GPU_ID -nt $NT \
     -chunks 5 25

bash rotate_CG_shifted.sh -f ../run/run/${NAME}25.gro -p ../run/${NAME}25_rotate_topol.top \
     -t ../run/run/${NAME}25.cpt --index ../run/${NAME}.ndx -deffnm ${NAME}45 -gpu_id $GPU_ID -nt $NT \
     -chunks 25 45
