#!/bin/bash

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

cat << EOF > annealing.mdp
title                   = Annealing with decreasing restraints
; Run parameters
integrator              = md                  ; leap-frog integrator
nsteps                  = 17000000              ; 34 ns total for annealing
dt                      = 0.002               ; 2 fs
; Output control
nstxout                 = 0                 ; save coordinates every 1.0 ps
nstvout                 = 0                 ; save velocities every 1.0 ps
nstenergy               = 500                 ; save energies every 10 ps
nstlog                  = 0                 ; update log file every 1.0 ps
compressed-x-grps       = non-Water ; Dont save water coordinates
nstxout-compressed      = 5000      ; save compressed coordinates avery 10 ps

; Bond parameters
continuation            = yes                 ; continuing from NVT
constraint_algorithm    = lincs               ; holonomic constraints 
constraints             = h-bonds             ; bonds involving H are constrained
lincs_iter              = 1                   ; accuracy of LINCS
lincs_order             = 4                   ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 40        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME                 ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4                   ; cubic interpolation
fourierspacing          = 0.16                ; grid spacing for FFT
; Temperature coupling with annealing
tcoupl                  = V-rescale           ; modified Berendsen thermostat
tc-grps                 = Other Water_and_ions ; two coupling groups
tau_t                   = 0.1   0.1           ; time constant, in ps
ref_t                   = 300   300           ; refernce temp
annealing               = single single       ; perform annealing for both groups
annealing-npoints       = 4 4                 ; 5 annealing points for both groups
annealing-time          = 0 1000 4000 5000 0 1000 4000 5000; time points for each group [ps]
annealing-temp          = 300 450 450 300 300 450 450 300 ; temp at each time point [K]
; Pressure coupling - starting with weak coupling
pcoupl                  = Berendsen           ; Berendsen barostat - gentler for annealing
pcoupltype              = semiisotropic       ; uniform scaling of x-y box vectors, independent z
tau_p                   = 1.0                 ; time constant, in ps
ref_p                   = 1.0 1.0             ; reference pressure, in bar
compressibility         = 4.5e-5 4.5e-5       ; isothermal compressibility of water
refcoord_scaling        = com                 ; scale COM of reference coordinates
; Periodic boundary conditions
pbc                     = xyz                 ; 3-D PBC
; Velocity generation
gen_vel                 = no                  ; continue from previous velocities
EOF


gmx grompp -f annealing.mdp -c $FILE_PATH -r $FILE_PATH -p ${TOPOLOGY_PATH} -o ../tmp/annealing_${OUTPUT_NAME}.tpr \
  -maxwarn 1 > ../logs/annealing_${OUTPUT_NAME}_grompp.out

gmx mdrun -deffnm ../tmp/annealing_${OUTPUT_NAME} -nt $NCORES -gpu_id $GPU_ID -bonded gpu -update gpu

# Run NPT equilibration
gmx grompp -f npt_equilibration.mdp -c ../tmp/annealing_${OUTPUT_NAME}.gro -p ${TOPOLOGY_PATH}  \
    -o ../tmp/annealing_npt_${OUTPUT_NAME}.tpr -t ../tmp/annealing_${OUTPUT_NAME}.cpt -maxwarn 1
gmx mdrun -deffnm ../tmp/annealing_npt_${OUTPUT_NAME} -nt $NCORES -gpu_id $GPU_ID -bonded gpu -update gpu 
echo 13 | gmx editconf -f ../tmp/annealing_npt_${OUTPUT_NAME}.tpr -o ../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}.pdb -n ../system/${OUTPUT_NAME}/index.ndx -conect

mkdir -p ../system/${OUTPUT_NAME}/
mv ../tmp/annealing*${OUTPUT_NAME}* ../system/${OUTPUT_NAME}/

bash util_annealing_plot.sh ${OUTPUT_NAME}
