#!/bin/bash/
mkdir -p ../tmp/
mkdir -p ../logs/
mkdir -p ../structure_files/
rm step*
rm \#*
rm bck*
rm ../*/\#*

TOPOLOGY_PATHS=()
NUMBER_MOLECULES=()
TOTAL_MOLECULES=0
N_STEPS_SCALE=21
NCORES=16
PIN=0

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--fiber_name)
      shift #past argument
      FILE_PATH=$1
      x=${1%.pdb}
      NAME=("${x##*/}")
      shift #past value
      ;;
    -t|--topologies)
      shift #past argument
      while [[ $1 == *.top ]]; do
         TOPOLOGY_PATHS+=("$1")
         shift #past value
      done
      ;;
    -p|--topology)
      shift #past argument
      TOPOLOGY_PATH=$1
      shift #past value
      ;;
    -ns|--n_scale)
      shift #past argument
      N_STEPS_SCALE=$((10+$1))
      shift #past value
      ;;
    -n|--number_of_molecules)
      shift #past argument
      while [[ $1 =~ ^[0-9]+$ ]]; do
         NUMBER_MOLECULES+=("$1")
         TOTAL_MOLECULES=$(($TOTAL_MOLECULES+$1))
         shift #past value
      done
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

bash build_topology.sh -t ${TOPOLOGY_PATHS[@]} -n ${NUMBER_MOLECULES[@]} -ion NA -nc $(($TOTAL_MOLECULES*2)) -o $TOPOLOGY_PATH
cp  $TOPOLOGY_PATH ../system/${NAME}/${NAME}.top

gmx editconf -f $FILE_PATH -o ../tmp/distbox_${NAME}.gro -d 0.1
z=$(tail -n 1 ../tmp/distbox_${NAME}.gro | awk '{print $3}')
gmx editconf -f ../tmp/distbox_${NAME}.gro -o ../tmp/em_${NAME}.gro -box 100 100 $z

rm ../tmp/\#*
rm step*
cat > emin3.mdp << EOF
integrator      = steep     ; Algorithm (steep = steepest descent minimization)
emtol           = 1000.0    ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep          = 0.01      ; Energy step size
nsteps          = 50000     ; Maximum number of (minimization) steps to perform

nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
rlist           = 1.2       ; Cut-off for making neighbor list (short range forces)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2       ; Short-range electrostatic cut-off
rvdw            = 1.2       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions
EOF

for i in $(seq 0 10); do
    echo "Iteration $i"
    . scale.sh -s 0.95 -f ../tmp/em_${NAME}.gro -ax xy
    gmx grompp -f emin3.mdp -c ../tmp/em_${NAME}.gro -p $TOPOLOGY_PATH -o ../tmp/em_${NAME}.tpr -maxwarn 1 
    gmx mdrun -deffnm ../tmp/em_${NAME} -nt 16 -gpu_id $GPU_ID
    echo $(tail -n 1 ../tmp/em_${NAME}.gro)
    cp ../tmp/em_${NAME}.gro ../tmp/em_${NAME}_$i.gro
done

for i in $(seq 11 $N_STEPS_SCALE); do
    echo "Iteration $i"
    echo 0 | gmx traj -f ../tmp/em_${NAME}.gro -s ../tmp/em_${NAME}.tpr -oxt ../tmp/em_whole_${NAME}.gro -pbc > /dev/null 2>&1
    . scale.sh -s 0.95 -f ../tmp/em_whole_${NAME}.gro -ax z
    gmx grompp -f emin3.mdp -c ../tmp/em_whole_${NAME}.gro -p $TOPOLOGY_PATH -o ../tmp/em_${NAME}.tpr -maxwarn 1 > /dev/null 2>&1
    gmx mdrun -deffnm ../tmp/em_${NAME} -nt 16 -gpu_id $GPU_ID
    echo $(tail -n 1 ../tmp/em_${NAME}.gro)
    cp ../tmp/em_${NAME}.gro ../tmp/em_${NAME}_$i.gro
done

cat > npt_short.mdp << EOF
integrator              = md        ; leap-frog integrator
nsteps                  = 10     ; 1000 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0       ; save coordinates every 1.0 ps
nstvout                 = 0       ; save velocities every 1.0 ps
nstenergy               = 0       ; save energies every 1.0 ps
nstlog                  = 0       ; update log file every 1.0 ps
nstxout-compressed      = 0      ; save compressed coordinates every 1.0 ps
; Bond parameters
continuation            = no       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = System
tau_t                   = 0.5        ; time constant, in ps
ref_t                   = 300       ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Berendsen     ; Pressure coupling on in NPT
pcoupltype              = semiisotropic             ; uniform scaling of box vectors
tau_p                   = 0.5                   ; time constant, in ps
ref_p                   = 0 1.0                   ; reference pressure, in bar
compressibility         = 0 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 1-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 
EOF
gmx grompp -f npt_short.mdp -c ../tmp/em_${NAME}.gro -p $TOPOLOGY_PATH -o ../tmp/npt_short_${NAME}.tpr -maxwarn 1 > /dev/null 2>&1
gmx mdrun -deffnm ../tmp/npt_short_${NAME} -nt 16 -gpu_id $GPU_ID

mkdir -p ../system/${NAME}
cp ../tmp/npt_short_${NAME}.gro ../system/${NAME}/${NAME}_dense.gro
cp $TOPOLOGY_PATH ../system/${NAME}/${NAME}.top
