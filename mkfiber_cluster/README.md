# mkfiber_cluster

## Purpose
- Builds coarse-grained bundles of multiple fibres (hexagonal 3×N layers) to
  probe inter- and intrafibrillar properties: diffusion coefficients,
  radial densities, ion coordination, SOAP descriptors, etc.  Extends the
  single-fibre `mkfiber` workflow to crowded bundles and produces data used in
  the bundle section of the manuscript.
- Inputs are the CG fibre structures/force fields emitted by `mkforce_field`
  and `mkfiber`, while outputs include equilibrated bundle trajectories, density/diffusion CSVs, and figures saved under
  `mkfiber_cluster/{run,system,logs,SOAP,figures}`.

## Inputs and Outputs
- **Inputs:**
  - CG fibre reference structure `structure_files/fibre_structure.pdb`, CG
    topology `structure_files/CG_system.top`, and counter-ion PDBs.
- **Outputs:**
  - Bundle build/equilibration intermediates under `system/<name>/`
    (`compressed.gro`, `equilibrated.gro`, topologies/indices).
  - Production runs plus `SOAP/`, analysis, and derived density/diffusion tables in `figures/` or `data/`.

## How to Use

### Prerequisites
- Ensure CG fibre inputs exist in `structure_files/`.

### Build and run the 5-layer bundle
```bash
cd mkfiber_cluster/bin
bash run_all.sh 0 32
```
- Arguments: `<gpu_id> <nt>`.  The script:
  1. Calls `packmol_structure.sh` to arrange 5 layers of fibres plus ions.
  2. Uses `build_system.sh` and `equilibriate.sh` to compress and solvate the
     bundle (results written to `system/5_layers/`).
  3. Runs a production simulation (`run/run/5_layers`) and automatically invokes
     `compute_diffusion.sh` to generate diffusion metrics for the bundle.

### Analyse diffusion and density
```bash
cd mkfiber_cluster/bin
python compute_diffusion.py ../run/run/5_layers.xtc ../run/5_layers.ndx 5_layers
python compute_density.py ../run/run/5_layers.xtc ../run/5_layers.ndx 5_layers
```
- Scripts expect the cleaned trajectory and index file, generate CSVs under
  `../figures/` or `../data/`, and create helper plots/notebooks
  (`analysis.ipynb`, `compute_density.ipynb`, `soapify.ipynb`) for detailed
  interpretation.

## Known Issues
- Building fresh bundles via `packmol_structure.sh` is fragile; many packings
  collapse during compression or equilibration.  Stable starting points are kept
  in `system/<name>/equilibrated.gro` and should be reused when possible.
- Scripts assume a fixed 5-layer, 3-fibre-per-layer geometry (total fibres and
  ion counts are hard-coded).  Adapting to other lattice sizes requires editing
  the Packmol input generation and topology duplication logic manually.
