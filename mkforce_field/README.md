
# mkforce_field

## Purpose
- Generates all-atom force-field parameters (RESP charges, `.gro/.top/.pdb`
  files) for molecular motors so that every other module in
  `molecular_muscle` can assemble fibres or coarse-grained models from a
  consistent parameter set.
- Inputs come from `structure_files/mol2/` (fragment `.mol2` files plus
  manually curated `main_chain/*.chain` templates).  Outputs are copied to
  `structure_files/{gro,pdb,top}` and used by `mkfiber`,
  `mkfiber_cluster`, and `coarse_grain_*`.

## Inputs and Outputs
- **Inputs:**
  - Fragment geometries in `structure_files/mol2/` (`head.mol2`,
    `core<SUBTYPE>-{R,S}.mol2`, `right_leg<SUBTYPE>.mol2`, etc.).
  - Chain layout templates in `mkforce_field/main_chain/*.chain` describing how
    residues join (extend or duplicate these files for new layouts).
  - Charges supplied via `run_molecule.sh` flags (`--head_charge`, etc.); net
    charge may be non-zero if the downstream simulations require it.
- **Outputs:**
  - GROMACS-ready `.gro/.top/.pdb` files for each enantiomer placed in
    `structure_files/{gro,top,pdb}`.
  - Intermediate `.mol2`, `.prepi`, `.frcmod`, `.prmtop`, `.inpcrd` stored in
    `mkforce_field/output/` and `mkforce_field/tmp/` for debugging or reuse.
  - Amber/ParmEd logs in `mkforce_field/logs/`.

## How to Use

### Prerequisites
- AmberTools (`antechamber`, `prepgen`, `parmchk2`, `tleap`) on `PATH`.
- ParmEd (`parmed`) and GROMACS (`gmx`) available, typically via the
  `molecular-motors` conda environment referenced inside `run_molecule.sh`.
- Populate `structure_files/mol2/` with the fragment inputs and ensure the
  matching `main_chain/*.chain` definitions exist (copy/modify the provided
  ones for new chemistries).

### Reproduce published results
```bash
cd mkforce_field/bin
bash run_molecule.sh --sub_type 2m --head_charge 0 --core_charge 0 --leg_charge -1
```
- Loops over `R`/`S` isomers, runs `mol2_to_prepi.sh` + `build_sequence.sh`,
  converts Amber output to `.gro/.top`, aligns with `gmx editconf`, then copies
  final files into `structure_files/{gro,pdb,top}`.  Inspect
  `mkforce_field/logs/` for Amber output if a step fails.

### Additional workflows
**Parametrise a single fragment (e.g., core for coarse-graining reference)**
```bash
cd mkforce_field/bin
bash run_small_molecule.sh --file_name ../../structure_files/mol2/core2m-R --name core2m-R
```
- Performs `antechamber → parmchk2 → tleap → ParmEd` on the specified `.mol2`
  and drops the resulting `.gro/.top` pair into `mkforce_field/output/`.

**Custom fragment layout**
```bash
cd mkforce_field/bin
bash mol2_to_prepi.sh --file_name ../../structure_files/mol2/head.mol2 \
                      --name HEAD --main_chain ../main_chain/head.chain \
                      --residue_name HED --net_charge 0
bash build_sequence.sh HEAD CORE6-R RLEG6 \
  --sequence "HED COR RLG" --output motor6-R
```
- Example of manually parametrising fragments before sequencing them; edit the
  `--sequence` string to include additional residues or alternative layouts.

## Known Issues
- If new fragment naming schemes are introduced, remember to update the
  `main_chain/*.chain` templates; missing/incorrect dummy atoms manifest as
  `prepgen` errors.
