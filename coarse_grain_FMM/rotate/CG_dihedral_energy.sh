#!/bin/bash
NCORES=$1
GPU_ID=$2
mkdir -p CG-metad
cp ../match-dist/v4/CG_system.top CG-metad/CG_system.top
cp ../match-dist/v4/CG_motor2m-R_v4.itp CG-metad/CG_motor2m-R_v4.itp
cp ../match-dist/v4/run/run.gro CG-metad/CG_motor.gro
cp ../match-dist/v4/index.ndx CG-metad/index.ndx
cd CG-metad

#!/bin/bash
cat > plumed.dat << EOF
WHOLEMOLECULES ENTITY0=1-28

phi:  TORSION ATOMS=5,1,8,13
psi:  TORSION ATOMS=5,7,8,11

metad: METAD ARG=phi PACE=10 HEIGHT=0.3 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 
    GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 TEMP=300
    ...

PRINT STRIDE=10 ARG=* FILE=COLVAR
EOF

cat > run.mdp << EOF
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000000 ; 500 000 000 * 0.002fs = 100 ns (1 mus)
dt                      = 0.02     ; 20 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 500      ; save energies every 10.0 ps
nstlog                  = 500      ; update log file every 10.0 ps
nstxout-compressed      = 500      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = MOTOR    ; save the whole system

cutoff-scheme            = Verlet
nstlist                  = 30
ns_type                  = grid
pbc                      = xyz
verlet-buffer-tolerance  = 0.005

; Electrostatics
coulombtype              = Reaction-Field ; I DO NOT KNOW HOW THIS WORKS
coulomb-modifier         = Potential-shift-verlet
rcoulomb                 = 1.1
epsilon_r                = 15	; 2.5 (with polarizable water)
epsilon_rf               = 80
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = MOTOR SOL_ION   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 12.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 3e-4                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com

gen_vel                  = no
gen_temp                 = 300

constraints              = none 
constraint_algorithm     = Lincs
EOF

gmx grompp -p CG_system.top -c CG_motor.gro -n index.ndx -f run.mdp -o run.tpr 
gmx mdrun -v -deffnm run -plumed plumed.dat -nb cpu -nt $NCORES -gpu_id $GPU_ID

rm COLVAR-reweight phi-reweight.weight
cat << EOF > plumed_reweight.dat
RESTART
WHOLEMOLECULES ENTITY0=1-28

phi:  TORSION ATOMS=5,1,8,13
psi:  TORSION ATOMS=5,7,8,11

metad: METAD ARG=phi PACE=100000000 HEIGHT=0.3 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 
    GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 TEMP=300
    ...

PRINT STRIDE=1 ARG=phi,metad.bias FILE=COLVAR-reweight
EOF
plumed driver --mf_xtc run.xtc --plumed plumed_reweight.dat

bmax=`awk 'BEGIN{max=0.}{if($1!="#!" && $3>max)max=$3}END{print max}' COLVAR-reweight`

# print phi values and weights
awk '{if($1!="#!") print $2,exp(($3-bmax)/kbt)}' kbt=2.494339 bmax=$bmax COLVAR-reweight > phi-reweight.weight


rm \#*
rm */\#*
rm bck*
cd ..
