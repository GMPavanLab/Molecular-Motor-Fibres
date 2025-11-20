# coarse_grain_fibre

## Purpose
- Maps the all-atom (AA) fibre trajectories from `mkfiber` onto a coarse-grained
  (CG) representation, simulates the CG fibre in solvent, and compares structural
  observables (radius, SASA, Voronota volumes, g(r), diffusion) between the AA and CG
  ensembles.  
- Uses AA reference runs from `coarse_grain_fibre/AA-reference/` and produces CG
  trajectories in `CG-reference/`, with cleaned data/plots stored in `trajectories/`
  and `logs/`.

## Inputs and Outputs
- **Inputs:**
  - AA fibre trajectories/topologies (`AA-reference/AA-run*.xtc/.tpr/.top`) and
    the bead mapping file `structure_files/mapping.ndx`.
  - CG force-field assets (`structure_files/CG_system*.top`, water box gro files).
- **Outputs:**
  - Mapped CG trajectories (`mapp-traj/CG-fibre<name>.gro/.xtc`) plus hacked AA
    topologies/itps for debugging.
  - CG production runs (`CG-reference/run/run<name>.*`, `eq*/`, `em/`) and the
    corresponding indices/energy logs.
  - Cleaned trajectories and analysis CSVs in `trajectories/<name>/` and plots in
    `logs/<analysis>/`.

## How to Use

### Prerequisites
- Ensure the AA reference simulations in `AA-reference/` are complete and that
  `structure_files/mapping.ndx` reflects the desired bead mapping (currently 152
  AA atoms → 28 beads across 100 motors).

### Pipeline to compare AA vs CG
1. **Map the AA fibre to CG beads**
   ```bash
   cd coarse_grain_fibre/bin
   bash mapp-AA.sh 2m
   ```
   Outputs `mapp-traj/CG-fibre2m.gro/.xtc`, cleaned topologies, and generates a
   CG `.tpr` ready for simulation.

2. **Run the CG solvent simulation**
   ```bash
   cd coarse_grain_fibre/bin
   bash run_CG.sh 2m 0
   ```
   (arguments: subtype/name and GPU id).  The script solvates the mapped fibre,
   removes overlapping water, adds ions, stitches index groups, carries out EM +
   five equilibration stages, then a production run.  Results live in
   `CG-reference/`.

3. **Clean trajectories for analysis**
   ```bash
   cd coarse_grain_fibre/bin
   bash clean_trajectory.sh \
       -s ../CG-reference/run/run2m.gro \
       -f ../CG-reference/run/run2m.xtc \
       -p ../CG-reference/CG_system2m.top \
       -t ../CG-reference/run/run2m.tpr \
       -o cg_2m
   ```
   Creates `trajectories/cg_2m/cg_2m.{pdb,xtc,ndx}` analogous to the AA workflow.

4. **Run analyses (radius, diffusion, SASA, Voronota)**
   ```bash
   cd coarse_grain_fibre/bin
   python measure_parameters.py ../trajectories/cg_2m/cg_2m.pdb \
          ../trajectories/cg_2m/cg_2m.xtc cg_2m ALL
   bash compute_distr.sh
   bash compute_sasa.sh
   python compute_diffusion.py
   ```
   Logs and figures populate `../logs/` and `../SASA/`.

### Optional notebooks
- `analysis.ipynb`, `analysis_measurements.ipynb`, and `rotate/analyse_rotation.ipynb`
  load the CSVs produced above to reproduce the manuscript plots.

## Known Issues
- The mapping and CG simulation scripts assume a fixed fibre size: 100 AA motors
  mapped to 200 CG molecules with 28 beads each, split into iso-R/iso-S blocks.
  Different stoichiometries require manual edits to `mapp-AA.sh`, `run_CG.sh`,
  and the generated index selections.