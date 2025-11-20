# bend_fibre

## Purpose
- Builds CG fibre bundles with imposed rotations to probe how a longer,
  partially rotated fibre relaxes when embedded between straight neighbours.
  Automation here mirrors `mkfiber_cluster` but with a fixed 3-fibre
  arrangement where one fibre is rotated.
- Inputs reuse CG structures/force-fields from `structure_files/` and outputs
  include the compressed/equilibrated bundles (`system/<name>/`), and rotated
  production runs (`run/run/<name>*.{gro,xtc}`).

## Inputs and Outputs
- **Inputs:**
  - `structure_files/fibre_structure.pdb` or `.gro` plus `structure_files/CG_system.top`.
  - Ion PDB (`structure_files/CA.pdb`), solvent box (`structure_files/box_CG_W_eq.gro`)..
- **Outputs:**
  - Build intermediates and topologies in `system/<name>/` (`compressed.gro`,
    `equilibrated.gro`, `*_topol*.top`).
  - Rotated production trajectories in `run/run/<name>*`.

## How to Use

### Prerequisites
- Same environment as `mkfiber_cluster`
- Ensure CG reference files and solvent boxes exist in `structure_files/`.

### Build and equilibrate the bent bundle (shifted workflow)
```bash
cd bend_fibre/bin
bash run_all.sh 0 32 bend_demo
```
- Arguments: `<gpu_id> <nt> <name>`.  The script:
  1. Invokes `packmol_structure_shifted.sh` to pack three fibres whose centers
     are shifted along z so one fibre can be rotated independently.
  2. Runs compression (`build_system.sh`) and the shifted solvation/equilibration
     pipeline (`equilibriate_shifted.sh`) to populate `system/<name>/`.
  3. Applies staged PLUMED rotations with `rotate_CG_shifted.sh` (chunks 0–5,
     5–25, 25–45 by default), producing intermediate trajectories under
     `run/run/<name>*`.

- The original symmetric workflow (`packmol_structure.sh`,
  `equilibriate.sh`, `rotate_CG.sh`) remains available for experiments that do
  not need shifted fibres or to analyse self-assembly of edges.

### Apply staged rotations manually
```bash
cd bend_fibre/bin
bash rotate_CG_shifted.sh \
    -f ../system/bend_demo/equilibrated.gro \
    -p ../system/bend_demo_topol_sol.top \
    -t ../system/bend_demo/equilibrated.cpt \
    --index ../run/bend_demo.ndx \
    -deffnm bend_demo_rot \
    -gpu_id 0 -nt 32 \
    -chunks 0 50
```
- Reuses the equilibrated bundle and rotates a selected number of motors
  (`-chunks`). Note a chunk is a number of randomly distributed motors on the fibres, not a segment. Chain multiple calls (with different checkpoints and chunk
  ranges) to replicate the staged bending used in `run_all.sh`.

## Known Issues
- Rotation remains highly sensitive and very small chunks or several random seeds may be needed to achieve a viable rotation scheme. 
- Geometry parameters (3 fibres, specific ion counts) and rotation chunks are
  hard-coded; changing fibre counts or segment lengths requires editing
  Packmol/topology duplication logic and PLUMED loops manually.
- The rotation workflow assumes PLUMED seed values and chunk indices that match
  the original setup; new configurations may also need recalibrated `plumed.dat`
  definitions to avoid instabilities.
