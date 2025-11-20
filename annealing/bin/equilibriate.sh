#!/bin/bash
#Duplicates the length of the fibre, solvates and equillibrates a large fiber system
NUMBER_OF_DUPLICATES=2
NCORES=16
PIN=0
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
    -n|--number_of_duplicates)
      shift #past argument
      NUMBER_OF_DUPLICATES=$1
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
    -nsteps|--number_steps)
      NSTEPS=$2
      shift;shift
      ;;
    -pin)
      PIN=$2
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

#Duplicate system
mkdir -p ../system/
mkdir -p ../system/${OUTPUT_NAME}/
cp ${TOPOLOGY_PATH} ../tmp/${OUTPUT_NAME}_tmp.top
TOPOLOGY_PATH=../tmp/${OUTPUT_NAME}_tmp.top

gmx editconf -f $FILE_PATH -o ../tmp/distbox.gro -d 0
z=$(tail -n 1 $FILE_PATH | awk '{print $3}')
gmx editconf -f $FILE_PATH -o ../tmp/periodic.gro -box 10 10 $z

gmx genconf -f ../tmp/periodic.gro -o ../tmp/${OUTPUT_NAME}_box.gro -nbox 1 1 $NUMBER_OF_DUPLICATES
gmx genconf -f ../tmp/periodic.gro -o ../system/${OUTPUT_NAME}/${OUTPUT_NAME}.pdb -nbox 1 1 $NUMBER_OF_DUPLICATES 

sed -n '/\[ molecules \]/,$p' ${TOPOLOGY_PATH} | sed 1,2d >> ${TOPOLOGY_PATH} #Duplicate molecule numbers in top file

#Solvate system and remove water inside fiber
gmx solvate -cp ../tmp/${OUTPUT_NAME}_box.gro -cs spc216.gro -p ${TOPOLOGY_PATH} -o ../tmp/${OUTPUT_NAME}_solv.gro
python3 remove_water_inside_fiber.py ../tmp/${OUTPUT_NAME}_solv.gro ../tmp/${OUTPUT_NAME}_rmwater.gro ${TOPOLOGY_PATH} 3

#Run constrained energy minimisation
cat << EOF > "emin2.mdp"
title                   = Second energy minimization
define                  = -DPOSRES_FIBRE ; position restrain the heads of the fiber

integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF

gmx grompp -f emin2.mdp -c ../tmp/${OUTPUT_NAME}_rmwater.gro -r ../tmp/${OUTPUT_NAME}_rmwater.gro -p ${TOPOLOGY_PATH} -o ../tmp/em2_${OUTPUT_NAME}.tpr  -maxwarn 1 > ../logs/em2_${OUTPUT_NAME}_grompp.out
gmx mdrun -deffnm ../tmp/em2_${OUTPUT_NAME} -nt $NCORES -gpu_id $GPU_ID > ../logs/em2_${OUTPUT_NAME}_run.out
printf '11\n0\n' | gmx energy -f ../tmp/em2_${OUTPUT_NAME}.edr -o ../output/em2_${OUTPUT_NAME}_potential.xvg

#Run constrained NVT equilibration
cat << EOF > "nvt_equilibration1.mdp"
title                   = Constrained Fiber NVT equilibration 
define                  = -DPOSRES_FIBRE ; position restrain the fiber

; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0       ; save coordinates every 1.0 ps
nstvout                 = 0       ; save velocities every 1.0 ps
nstenergy               = 0       ; save energies every 1.0 ps
nstlog                  = 0       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Other Water_and_ions   ; coupling fiber and water+ions
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
EOF

gmx grompp -f nvt_equilibration1.mdp -c ../tmp/em2_${OUTPUT_NAME}.gro -r ../tmp/em2_${OUTPUT_NAME}.gro -p ${TOPOLOGY_PATH} -o ../tmp/nvt_${OUTPUT_NAME}.tpr  -maxwarn 1 > ../logs/em2_${OUTPUT_NAME}_grompp.out
gmx mdrun -deffnm ../tmp/nvt_${OUTPUT_NAME} -nt $NCORES -gpu_id $GPU_ID > ../logs/nvt_${OUTPUT_NAME}_run.out

#Run unconstrained NVT equilibration
#Run constrained NVT equilibration
cat << EOF > "nvt_equilibration2.mdp"
title                   = Constrained Fiber NVT equilibration 

; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0       ; save coordinates every 1.0 ps
nstvout                 = 0       ; save velocities every 1.0 ps
nstenergy               = 0       ; save energies every 1.0 ps
nstlog                  = 0       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 30        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Other Water_and_ions   ; coupling fiber and water+ions
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; Velocity generation is off
EOF

gmx grompp -f nvt_equilibration2.mdp -c ../tmp/nvt_${OUTPUT_NAME}.gro -p ${TOPOLOGY_PATH} -o ../tmp/nvt_${OUTPUT_NAME}2.tpr \
  -maxwarn 1 > ../logs/em2_${OUTPUT_NAME}_grompp.out

gmx mdrun -deffnm ../tmp/nvt_${OUTPUT_NAME}2 -nt $NCORES -gpu_id $GPU_ID > ../logs/nvt_${OUTPUT_NAME}2_run.out

#Run NPT equilibration
cat << EOF > "npt_equilibration.mdp"
title                   = Unconstrained Fibril NPT Equilibration
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 1000000     ; 2 ns
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0       ; save coordinates every 1.0 ps
nstvout                 = 0       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 0       ; update log file every 1.0 ps
compressed-x-grps       = non-Water ; Dont save water coordinates
nstxout-compressed      = 500000      ; save compressed coordinates 

; Bond parameters
continuation            = yes       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 40        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Other Water_and_ions   ; two coupling groups - more accurate
tau_t                   = 0.5     0.5           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Berendsen     ; Pressure coupling on in NPT
pcoupltype              = semiisotropic             ; uniform scaling of box vectors
tau_p                   = 1.0                   ; time constant, in ps
ref_p                   = 1.0 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 
EOF

gmx grompp -f npt_equilibration.mdp -c ../tmp/nvt_${OUTPUT_NAME}2.gro -t ../tmp/nvt_${OUTPUT_NAME}2.cpt -p ${TOPOLOGY_PATH} -o ../tmp/npt_long_${OUTPUT_NAME}.tpr \
  -maxwarn 1 > ../logs/npt_long_${OUTPUT_NAME}_grompp.out
cp ../tmp/npt_long_${OUTPUT_NAME}.tpr ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}.tpr
cp ${TOPOLOGY_PATH} ../system/${OUTPUT_NAME}.top

gmx mdrun -deffnm ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME} -nt $NCORES -gpu_id $GPU_ID -bonded gpu -update gpu
echo q | gmx make_ndx -f ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}.tpr -o ../system/${OUTPUT_NAME}/index.ndx
echo 13 | gmx editconf -f ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}.tpr -o ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}.pdb -n ../system/${OUTPUT_NAME}/index.ndx -conect
printf '21\n0\n' | gmx energy -f ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}.edr -o ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}_z.xvg
printf '23\n0\n' | gmx energy -f ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}.edr -o ../system/${OUTPUT_NAME}/npt_${OUTPUT_NAME}_density.xvg

rm ../tmp/*#*#*