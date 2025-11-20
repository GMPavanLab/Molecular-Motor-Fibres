# characterise_fibre

## Purpose
- Cleans and analyses fibre trajectories produced by `mkfiber` (or anywhere else).  The scripts
  here unwrap, center, and reindex the raw GROMACS outputs before extracting
  structural metrics (radius, SASA, Voronota volumes/areas, torsions, etc.) that
  the plotting notebooks and downstream coarse-grain validation workflows use.
- Inputs primarily come from `mkfiber/system/<name>/` (NPT `.pdb/.xtc/.top/.tpr`)
  and the cleaned data are stored under `characterise_fibre/trajectories/`.
  Logs/plots feed back into this module and are referenced by notebooks in the
  same folder.

## Inputs and Outputs
- **Inputs:**
  - Raw simulation files from `mkfiber/system/<name>/` (or any compatible
    GROMACS run): `npt_<name>.pdb`, `.xtc`, `.top`, `.tpr`.
  - `structure_files/vdwradii_AA.dat` for SASA calculations, Voronota binaries,
    and the motor topology information (assumes 152-atom motors with named
    subdomains).
- **Outputs:**
  - Cleaned trajectories in `characterise_fibre/trajectories/<name>/<name>.{pdb,xtc,top,ndx}`.
  - Metric CSVs/plots in `characterise_fibre/logs/<analysis>/`.
  - Optional aggregate tables/notebook outputs saved alongside the notebooks
    (`analysis.ipynb`, `soapify.ipynb`, etc.).

## How to Use

### Prerequisites
- Cleaned trajectories under `characterise_fibre/trajectories/`.  If starting
  from mkfiber raw outputs, run `bin/clean_trajectory.sh` first.

### Reproduce published results
```bash
cd characterise_fibre/bin
bash run_all.sh
```
- Iterates over every directory in `../trajectories/` containing matching
  `<name>.pdb/.xtc`, then invokes `measure_parameters.py ... ALL` to run every
  available analysis method for each system.  Logs and CSVs are emitted in
  `../logs/`.

### Additional workflows
**Clean a single mkfiber trajectory**
```bash
cd characterise_fibre/bin
bash clean_trajectory.sh \
    -s ../../mkfiber/system/frac_75/npt_frac_75.pdb \
    -f ../../mkfiber/system/frac_75/npt_frac_75.xtc \
    -p ../../mkfiber/system/frac_75/frac_75.top \
    -t ../../mkfiber/system/frac_75/npt_frac_75.tpr \
    -o frac_75
```
- Rebuilds whole molecules, recenters the fibre, renames residues (LLG/RLG/HED),
  and writes cleaned files into `../trajectories/frac_75/`.

**Run selective analyses**
```bash
cd characterise_fibre/bin
python measure_parameters.py ../trajectories/frac_75/frac_75.pdb \
       ../trajectories/frac_75/frac_75.xtc frac_75 radius \
       --args resnames=RLG,LLG
python measure_parameters.py ../trajectories/frac_75/frac_75.pdb \
       ../trajectories/frac_75/frac_75.xtc frac_75 voronota \
       --args resids=1,2 resnames=COR
```
- Useful when debugging specific metrics or tuning residue/atom selections.

## Known Issues
- `clean_trajectory.sh` and `make_ndx.sh` assume each motor contains 152 atoms
  laid out like the mkfiber outputs; differing topologies require editing
  `NATOMS` and the subdomain ranges.
- Voronota-based analyses can be slow for long trajectories (processing every
  10th frame mitigates this, but runtimes remain on the order of hours).
