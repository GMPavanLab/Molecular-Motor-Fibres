#!/bin/bash
#Generates a geometry with motors in cylinder. 4 motors in height, $N motors in circonference. 
#Requires packmol and gromacs

FIBRE=$1
OUTPUT_FILE_NAME=$2
IONTYPE=$3
TOPOLOGY=$4
NLAYER=$5
GPU_ID=$6
OUTERR=40 #DERIVED ANALYTICALLY FROM d_opt monomer=11nm
INNERR=27

#Quick equilibration to avoid issues when removing counter ions
gmx grompp -f mdp/em.mdp -p $TOPOLOGY -c $FIBRE -o ../tmp/em.tpr -maxwarn 1 
gmx mdrun -deffnm ../tmp/em -v -nt 16 -pin on -pinoffset 32 -gpu_id $GPU_ID

cat > "../tmp/qeq.mdp" << EOF
integrator              = md        ; leap-frog integrator
nsteps                  = 10000     ; 0.02 * 500 = 10000 ps
dt                      = 0.02     ; 20 fs

cutoff-scheme            = Verlet
nstlist                  = 30
ns_type                  = grid
pbc                      = xyz
verlet-buffer-tolerance  = 0.005

coulombtype              = Reaction-Field ; I DO NOT KNOW HOW THIS WORKS
coulomb-modifier         = Potential-shift-verlet
rcoulomb                 = 1.1
epsilon_r                = 15	; 2.5 (with polarizable water)
epsilon_rf               = 80
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = System   ; two coupling groups - more accurate
tau_t                   = 0.1         ; time constant, in ps
ref_t                   = 100        ; reference temperature, one for each group, in K

pcoupl                  = Berendsen     ; Pressure coupling on in NPT
pcoupltype              = semiisotropic             ; uniform scaling of box vectors
tau_p                   = 4.0                   ; time constant, in ps
ref_p                   = 1.0 1.0                  ; reference pressure, in bar
compressibility         = 3e-4 3e-4               ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com

gen_vel                  = no
gen_temp                 = 100
gen_seed                 = 42

constraints              = none 
constraint_algorithm     = Lincs
EOF

gmx grompp -f ../tmp/qeq.mdp -p $TOPOLOGY -c ../tmp/em.gro -o ../tmp/qeq.tpr -maxwarn 1
gmx mdrun -deffnm ../tmp/qeq -v -nt 16 -pin on -pinoffset 32 -gpu_id $GPU_ID

gmx editconf -f ../tmp/qeq.gro -o ../tmp/fibre_vacuum.pdb -c
FIBRE=../tmp/fibre_vacuum.pdb
boxz=$(tail -n 1 ../tmp/qeq.gro | awk '{print $3}')
BOXLEN=$(echo "scale=2; 10*($boxz / 2)" | bc)
DIST=70
echo "BOXLEN: $BOXLEN"

mkdir -p ../tmp

cat > 'packmol.inp' << EOF

# All the atoms from diferent molecules will be at least 11 Angstrons
# apart at the solution

tolerance 1
# The output file name

output ../tmp/packmol.pdb
filetype pdb
writeout -1
EOF

for (( i=0; i<$NLAYER; i++)); do
cat >> 'packmol.inp' << EOF

# 3 fibres in hexagonal packing.
structure $FIBRE
  number 1
  center
  fixed $(echo "scale=2; ($DIST / 2 )*$(($i % 2))" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) 0. 0. 0. 0.
end structure 
structure $FIBRE
  number 1
  center
  fixed $(echo "scale=2; ($DIST / 2 )*$(($i % 2)) + $DIST" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) 0. 0. 0. 0.
end structure 
structure $FIBRE
  number 1
  center
  fixed $(echo "scale=2; ($DIST / 2 )*$(($i % 2))-$DIST" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) 0. 0. 0. 0.
end structure 
EOF
done
for (( i=0; i<$NLAYER; i++)); do
cat >> 'packmol.inp' << EOF
structure $IONTYPE
  number 220
  inside box -500 -500 -$BOXLEN 500 500 $BOXLEN
  outside cylinder $(echo "scale=2; ($DIST / 2 )*$(($i % 2))-$DIST" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) -200. 0. 0. 1. $INNERR 400
  inside cylinder $(echo "scale=2; ($DIST / 2 )*$(($i % 2))-$DIST" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) -200. 0. 0. 1. $OUTERR 400
end structure
structure $IONTYPE
  number 220
  inside box -500 -500 -$BOXLEN 500 500 $BOXLEN
  outside cylinder $(echo "scale=2; ($DIST / 2 )*$(($i % 2))" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) -200. 0. 0. 1. $INNERR 400
  inside cylinder $(echo "scale=2; ($DIST / 2 )*$(($i % 2))" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) -200. 0. 0. 1. $OUTERR 400
end structure
structure $IONTYPE
  number 220
  inside box -500 -500 -$BOXLEN 500 500 $BOXLEN
  outside cylinder $(echo "scale=2; ($DIST / 2 )*$(($i % 2)) + $DIST" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) -200. 0. 0. 1. $INNERR 400
  inside cylinder $(echo "scale=2; ($DIST / 2 )*$(($i % 2)) + $DIST" | bc) -$(echo "scale=2; $DIST*sqrt(1- 1/4)*$i" | bc) -200. 0. 0. 1. $OUTERR 400
end structure
EOF
done

#Run packmol on file
packmol < packmol.inp 

gmx editconf -f ../tmp/packmol.pdb -o ../tmp/tmp.gro -rotate 90 0 0
gmx editconf -f ../tmp/tmp.gro -o ../tmp/distbox.gro -d 0.1 -c
xlen=$(tail -n 1 ../tmp/distbox.gro | awk '{print $1}')
ylen=$(tail -n 1 ../tmp/distbox.gro | awk '{print $2}')
zlen=$(tail -n 1 ../tmp/distbox.gro | awk '{print $3}')
gmx editconf -f ../tmp/distbox.gro -o $OUTPUT_FILE_NAME -box $(echo "$ylen + 2" | bc -l) $boxz $(echo "6*$NLAYER*2" | bc -l)
cp $TOPOLOGY ../tmp/topol.top
for (( i=1; i<$(($NLAYER*3)); i++)); do
sed -n '/\[ molecules \]/,$p' ../structure_files/CG_system.top | sed 1,2d >> ../tmp/topol.top #Duplicate molecule numbers in top file
done
printf "CA $((220*3*$NLAYER))\n" >> ../tmp/topol.top #Add ions to top file
cp ../tmp/topol.top $TOPOLOGY