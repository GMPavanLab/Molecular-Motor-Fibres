# coarse_grain_MM2

## Purpose
- All-atom -> coarse-grained parametrisation workflow for the MM2
  motor core (and its actuated variant). Used to iteratively genenerate, test and tune new CG models.
- The subdir contains several drafts for models with bent stators, different mappings and parametrisations. Should not be used for production.
- Inputs come from `structure_files/` (AA structures, CG mappings, solvent
  boxes, index files), while outputs populate the subfolders below.
- To generate a new mapping, update the `structure_files/mapping.ndx`, then create a new .itp file.

## Inputs and Outputs
- **Inputs:**
  - AA reference structures (`structure_files/MM2*.gro/.top`), mapping files
    `structure_files/mapping_MM2*.ndx`, coarse-grained base topologies, and CG
    solvent boxes.
- **Outputs:**
  - `AA-reference/run/` – 1 µs AA trajectories for each motor variant.
  - `mapp-traj/CG-*.gro/.xtc/.tpr` – mapped CG trajectories derived from the AA
    references.
  - `match-dist/<version>/` – CG simulations used to fit bonds/angles/dihedrals,
    with distributions and comparison plots.
  - `SASA/AA` & `SASA/CG` – SASA traces/distributions for AA vs CG.
  - `check-parms/`, `rotate/`, `stacking/` – specialised analyses (parameter
    cross-checks, rotation tuning, dimer stacking PMFs).

## How to Use

### Prerequisites
- Populate `structure_files/` with AA structures (`MM2.gro/.top`,
  `MM2b-R.gro/.top`), CG solvent boxes, and mapping files
  (`mapping_MM2.ndx`, `mapping_MM2b-R.ndx`).

### Typical workflow
1. **Generate AA reference trajectories**
   ```bash
   cd coarse_grain_MM2/AA-reference
   bash run_all.sh 32 0 MM2
   bash run_all.sh 32 0 MM2b-R
   ```
   Produces solvated/em/nvt/npt and 1 µs production trajectories in
   `run/md_MM2*.xtc`.

2. **Map AA trajectories onto the CG beads**
   ```bash
   cd coarse_grain_MM2/mapp-traj
   bash run_all.sh MM2
   bash run_all.sh MM2b-R
   ```
   Outputs `mapped-MM2*.xtc`, `CG-MM2*.gro/.tpr` for use in CG fitting.

3. **Iteratively fit coarse-grained force fields**
   - Prepare the `[ bonds | angles | dihedrals ].ndx` you want to optimise in `structure_files/`.
   - Compute AA reference distributions:
     ```bash
     cd coarse_grain_MM2/match-dist
     . compute_distr.sh ../mapp-traj/mapped-MM2.xtc ../mapp-traj/CG-MM2.tpr reference
     ```
   - Create a topology candidate (`version*/CG_system_base.top`).
   - Run the CG simulation and extract distributions:
     ```bash
     bash run_version.sh --version v4 --molecule MM2 --gpu_id 0 --number_cores 32
     ```
   - Compare to AA distributions via `analysis.ipynb`, update topology, repeat.

4. **Measure SASA**
   ```bash
   cd coarse_grain_MM2/SASA
   bash compute_sasa.sh MM2 v4 CG
   bash compute_sasa.sh MM2 dummy AA
   ```
   Generates SASA time series for AA and CG references.

5. **Auxiliary analyses**
   - `check-parms/run_all.sh` & `analysis.ipynb`: verify AA parameters vs literature.
   - `rotate/`: PLUMED-driven tests for rotating the motor core.
   - `stacking/`: PMF calculations for dimerisation.

## Known Issues
- Pipeline assumes the legacy mapping that treats MM2b-R with the same bead
  layout as MM2 (only atom numbering differs).  Updating to a new mapping
  requires regenerating the `mapping_*.ndx` files and re-running the mapping/fit
  steps.
- Some subfolders (`rotate/`, `stacking/`) are marked “not up to date” and may
  require manual fixes before use.
- `run_version.sh` uses fixed solvent box dimensions and ion numbers; adjust
  manually if you deviate from the reference setup.
