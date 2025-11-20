#!/bin/bash
#Cant really remember why I made this script, but it seems to compute the dihedral energy for the AA system during rotation 
NCORES=$1
GPU_ID=$2
#AA dihedral
mkdir -p AA-metad
cp ../AA-reference/topol2m-R.top AA_system.top
cp ../AA-reference/run/md2m-R.gro AA-metad/AA_motor.gro
cp ../structure_files/index-R.ndx AA-metad/index.ndx
cd AA-metad
rm COLVAR phi.weight phi-reweight.weight COLVAR-reweight

#!/bin/bash
cat > plumed.dat << EOF
WHOLEMOLECULES ENTITY0=1-152

com5: COM ATOMS=117,118
com1: COM ATOMS=66,67,109
com7: COM ATOMS=110,111,116
com8: COM ATOMS=63,64
com11: COM ATOMS=39,40,51
com13: COM ATOMS=43,44,45,50

phi_CG: TORSION ATOMS=com5,com1,com8,com13
psi_CG: TORSION ATOMS=com5,com7,com8,com11
phi: TORSION ATOMS=63,64,65,110

metad: METAD ARG=phi PACE=100 HEIGHT=0.3 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 
    GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 TEMP=300
    ...

PRINT STRIDE=100 ARG=phi,phi_CG,psi_CG,metad.bias FILE=COLVAR
EOF

cat > run.mdp << EOF
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000000    ; 500 000 000 * 0.002fs = 1000 ns (1 mus)
dt                      = 0.002     ; 1 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = MOTOR    ; save the whole system
; Bond parameters
continuation            = no       ; Restarting after NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = MOTOR SOL_Ion   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 3e-4                ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correctionq
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 

EOF

gmx grompp -p ../AA_system.top -c AA_motor.gro -n index.ndx -f run.mdp -o run.tpr 
gmx mdrun -v -deffnm run -plumed plumed.dat -nt $NCORES -gpu_id $GPU_ID

rm COLVAR-reweight phi-reweight.weight
cat << EOF > plumed_reweight.dat
RESTART
WHOLEMOLECULES ENTITY0=1-152

com5: COM ATOMS=117,118
com1: COM ATOMS=66,67,109
com7: COM ATOMS=110,111,116
com8: COM ATOMS=63,64
com11: COM ATOMS=39,40,51
com13: COM ATOMS=43,44,45,50

phi_CG: TORSION ATOMS=com5,com1,com8,com13
psi_CG: TORSION ATOMS=com5,com7,com8,com11
phi: TORSION ATOMS=63,64,65,110

metad: METAD ARG=phi PACE=10000000 HEIGHT=0.3 SIGMA=0.05 ...
    FILE=HILLS BIASFACTOR=10 
    GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 TEMP=300
    ...
PRINT STRIDE=1 ARG=phi,phi_CG,psi_CG,metad.bias FILE=COLVAR-reweight
EOF

echo 1 | gmx editconf -f AA_motor.gro -o AA_motor.pdb -n index.ndx
plumed driver --mf_xtc run.xtc --plumed plumed_reweight.dat --pdb AA_motor.pdb

bmax=`awk 'BEGIN{max=0.}{if($1!="#!" && $5>max)max=$5}END{print max}' COLVAR-reweight`

# print phi values and weights
awk '{if($1!="#!") print $2,$3,$4,exp(($5-bmax)/kbt)}' kbt=2.494339 bmax=$bmax COLVAR-reweight > phi-reweight.weight



rm \#*
rm */\#*
rm bck*
cd ..
