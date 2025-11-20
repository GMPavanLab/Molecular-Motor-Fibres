# mkfiber

## Purpose
- Full pipeline to build, compress, solvate, and equilibrate all-atom fibres made from the
  motor `.pdb/.top` pairs created in `mkforce_field`, producing trajectories and
  dense structures that downstream analysis (`characterise_fibre`,
  `coarse_grain_fibre`, `mkfiber_cluster`) rely on.
- Reads structures from `structure_files/{pdb,top}` and pack-in parameters
  stored alongside the scripts in `mkfiber/bin`.  Outputs land in
  `mkfiber/system/<deffnm>/` (full trajectories + checkpoints) and
  `mkfiber/output/<deffnm>/` (dense `.gro/.top` pairs and equilibration logs).

## Inputs and Outputs
- **Inputs:**
  - Motor structures: `structure_files/pdb/motor<SUBTYPE>-{R,S}.pdb` (and `B`
    isomers when rotated motors are used).
  - Topologies: `structure_files/top/motor<SUBTYPE>-{R,S}.top`.
  - Simulation templates: Packmol input plus `.mdp` snippets in `mkfiber/bin`
    (scripts auto-write specific MDPs into `mkfiber/bin/mdp/`).
  - GPU/CPU resources accessible to GROMACS (`gmx`) and Packmol (`packmol`).
- **Outputs:**
  - Packed and compressed coordinates in `mkfiber/output/<name>/<name>_dense.gro`
    and `.top`, along with `.xvg` curves (`output/NPT/`) monitoring equilibration.
  - Fully solvated, doubled-length fibres plus trajectories/checkpoints in
    `mkfiber/system/<name>/` (NPT `.xtc`, `.cpt`, `.gro`, `.edr`).
  - Logs in `mkfiber/logs/` covering Packmol, GROMACS, and equilibration steps.

## How to Use

### Prerequisites
- Ensure `mkforce_field` has produced matching motor `.pdb/.top` pairs in
  `structure_files/`.

### Reproduce published results
```bash
cd mkfiber/bin
bash run_all.sh --subtype 2m --gpu_id 0
```
- Generates fibres across several rotation fractions (`run_frac_rot.sh`) and,
  with `frac=0`, repeats the pipeline using alternative compression scales
  (9/12/15) to emulate varying initial packings.  Each run invokes
  `make_fiber.sh`, which chains `packmol_structure.sh`, `compress_fibre.sh`, and
  `equilibriate_large_fiber.sh`.  Final metrics are extracted via
  `run_analysis.py`.

### Additional workflows
**Single rotation-fraction sweep**
```bash
cd mkfiber/bin
bash run_frac_rot.sh --subtype 2m --gpu_id 0 --frac 90 --deffnm frac_90 --n_scales 9
```
- Mixes straight and rotated motors according to `--frac`, executes the full
  make_fiber pipeline, and saves results under `mkfiber/output/frac_90`.

**Manual fibre assembly**
```bash
cd mkfiber/bin
. make_fiber.sh \
    --pdb_files ../../structure_files/pdb/motor2m-R.pdb ../../structure_files/pdb/motor2m-S.pdb \
    --top_files ../../structure_files/top/motor2m-R.top ../../structure_files/top/motor2m-S.top \
    --number_of_molecules 50 50 --deffnm demo --gpu_id 0
```
- Useful when experimenting with specific stoichiometries or ionisation states;
  call the helper scripts directly if you need to tune Packmol or compression.

## Known Issues
- `make_fiber.sh` assumes every motor carries a net charge of -2e when adding
  counter-ions during solvation; adjust the relevant commands manually if you
  introduce motors with different net charges.
- Some Packmol random seeds generate unstable initial fibres. Even with multiple
  compression scales, some systems may require manual re-packing or continued
  equilibration from saved `.tpr` checkpoints in `mkfiber/system/`. 
