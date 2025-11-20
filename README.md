# Molecular Motor Fibre Toolkit

Living muscles amplify nanometre-scale actin–myosin motion into centimetre-scale
work.  Artificial molecular muscles based on Feringa-type rotary motors follow
the same principle: amphiphilic motors self-assemble into worm-like nanofibres,
bundle into hydrated strings, and bend toward UV light when millions of motors
rotate in concert.  Small chemical tweaks, however, can disrupt this hierarchy,
and connecting molecular design to macroscopic actuation still demands a
multiscale modelling pipeline.  This repository delivers that pipeline: starting
from all-atom force-field derivation, through fibre assembly, coarse-grained
mapping, bundle simulations with synthetic rotations, and finally MDAnalysis-
based characterisation.

Each experimental phase lives in its own subdirectory with a README documenting
local tooling and automation scripts.  Combine them to trace how a molecular
motor’s chemistry propagates to fibre-scale mechanics and bundle response.

## Environments and Requirements
- `environment.yml` – Conda environment used for analysis notebooks and
  MDAnalysis-based scripts.
- `gromacs_info.txt` / `plumed_info.txt` – versions used for every production
  simulation (GROMACS 2021-modified with CUDA 11.6, PLUMED 2.7.1 built with MPI).
  The MDAnalysis tooling was compiled against PLUMED 2.9.2 and GROMACS 2024.5
  for post-processing; reproduce those builds if you need to re-run the
  notebooks.

## Project Layout
| Directory | Description | README |
| --- | --- | --- |
| `mkforce_field/` | Amber/ParmEd wrappers that derive AA force fields for new motor variants. | [mkforce_field/README](mkforce_field/README.md) |
| `mkfiber/` | Builds single fibres from AA motors, compresses, solvates, equilibrates. | [mkfiber/README](mkfiber/README.md) |
| `characterise_fibre/` | Cleans mkfiber trajectories and runs MDAnalysis/Voronota analysis. | [characterise_fibre/README](characterise_fibre/README.md) |
| `mkfiber_cluster/` | Packmol + CG bundle builder for 5-layer fibre packs and diffusion studies. | [mkfiber_cluster/README](mkfiber_cluster/README.md) |
| `bend_fibre/` | Generates shifted/rotated 3-fibre bundles to study bent motors. | [bend_fibre/README](bend_fibre/README.md) |
| `coarse_grain_fibre/` | AA vs CG fibre comparison, SASA, diffusion of bundled systems. | [coarse_grain_fibre/README](coarse_grain_fibre/README.md) |
| `coarse_grain_MM2/` | Full AA→CG pipeline for the MM2 motor core (mapping, fitting, SASA). | [coarse_grain_MM2/README](coarse_grain_MM2/README.md) |
| `coarse_grain_FMM/` | Functional-group analogue of the MM2 pipeline (MM2a motors). | [coarse_grain_FMM/README](coarse_grain_FMM/README.md) |
| `dimerisation PMF/` | PMF calculations for AA & CG motor dimers. | (see folder notes) |
| `annealing/` | Annealing protocol tests for AA fibres. | (experimental) |
| `structure_files/` | Shared AA/CG structures, topologies, mappings (see below). | — |

Each README documents the local `run_all` scripts, inputs, outputs, and known
issues.  Follow the sequence `mkforce_field → mkfiber → characterise_fibre`
before branching into the coarse-graining or bundle studies.

## Top-Level `run_all.sh`
```text
1. mkforce_field/bin/run_molecule.sh  # regenerates AA force fields
2. mkfiber/bin/run_all.sh             # builds and equilibrates fibres
3. characterise_fibre/bin/run_all.sh  # cleans & analyses the trajectories
```
Set `MOLTYPE` and `GPU_ID` inside `run_all.sh` to control which motor subtype
is built and which GPU to target.  Each internal script has additional flags for
fine-grained control if you need to deviate from the default pipeline.

## `structure_files/` Overview
- `gro/`, `pdb/`, `top/` – canonical AA structures and topologies for the motor
  fragments, used by mkforce_field and mkfiber.
- `mol2/` – source fragments for AmberTools (head/core/legs) in various
  subtypes.
- `number_elements_by_residue.py`, `sort_mol2_bonds.pl` – helper scripts for
  pre-processing input files.
- Mapping/CG assets (e.g., `mapping_MM2*.ndx`, CG solvent boxes) live alongside
  the directories listed above and are consumed by the coarse-graining modules.

## Contact
Questions or bug reports: **Mattias Karlsson** – `mattias.karlsson@example.com`
