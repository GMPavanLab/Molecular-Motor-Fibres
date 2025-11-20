# rotate_motors_in_fibre

## Purpose
- Applies PLUMED-MOVINGRESTRAINT to subsets of motors inside an
  all-atom fibre to study how enforced rotations affect structural stability
- Uses fibre initial conditions from `./structure_files/` (e.g.,
  `fibre_run.gro`, `fibre_structure.pdb`, `temp_groups.ndx`) and produces
  trajectories in `run/` plus COM tracks/COLVAR files in `positions/` and
  `COLVARS/`.

## Inputs and Outputs
- **Inputs:**
  - `structure_files/fibre_structure.pdb`, `fibre_run.gro`, `fibre_run.cpt`,
    `fibre.top`, and `temp_groups.ndx` defining the AA fibre with 100 motors.
  - PLUMED template `plumed.dat` (overwritten each run) and GROMACS MDP files
    (`mdrun_fix.mdp`).
  - User parameters passed to `run_AA.sh` (`NCORES`, `GPU_ID`, `KAPPA`,
    `tropic` label, `NROT` motors to rotate).
- **Outputs:**
  - Production trajectories/checkpoints in `run/mdrun<KAPPA>_<label>.{tpr,gro,xtc,edr}`.
  - PLUMED COLVAR files under `COLVARS/` and COM coordinate dumps in
    `positions/`.
  - Cleaned trajectories in `trajectories/` after running the
    `clean_trajectory.sh → measure_parameters.py` pipeline.

### Prerequisites
- Ensure the reference fibre files exist in `mdruns/rotate_motors_in_fibre/structure_files/`.

### Rotate a block of motors (article workflow)
```bash
cd mdruns/rotate_motors_in_fibre/bin
bash run_AA.sh 16 0 150 tropic_label 45
```
- Rotates `NROT` motors in both halves of the fibre using PLUMED moving
  restraints with stiffness `KAPPA`.  Outputs land in `../run/`,
  `../COLVARS/`, and `../positions/`.

### Clean and analyse a finished trajectory
```bash
cd mdruns/rotate_motors_in_fibre/bin
bash clean_trajectory.sh \
    -s ../run/mdrun150_tropic_label.pdb \
    -f ../run/mdrun150_tropic_label.xtc \
    -p ../run/topol.top \
    -t ../run/mdrun150tropic_label.tpr \
    -o tropic_label
python measure_parameters.py ../trajectories/tropic_label/tropic_label.pdb \
       ../trajectories/tropic_label/tropic_label.xtc tropic_label ALL
```
- Produces cleaned files in `../trajectories/tropic_label/` and runs the
  complete analysis suite (radius, torsions, Voronota, SASA, etc.).

## Known Issues
- The scripts assume exactly 100 motors with 152
  atoms each.  Systems built with other stoichiometries or atom counts require
  manual edits to `NMOLS`, `NATOMS`, and multiple atom index offsets inside
  `run_AA.sh`/`rotate_one.sh`.
- `clean_trajectory.sh` expects the same residue layout and indexing scheme as
  mkfiber outputs; custom systems may need updated `make_ndx.sh` logic before
  cleaning/analyzing.
