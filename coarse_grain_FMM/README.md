# coarse_grain_FMM

This directory mirrors the MM2 coarse-graining workflow but focuses on the
functional groups (MM2a).  Follow the detailed instructions in
`../coarse_grain_MM2/README.md` for the AA-reference generation, mapping,
match-dist fitting, and SASA analyses—the same tooling and assumptions apply.

Differences specific to the FMM folder:
- Only the functionalised cores (MM2a-R/S variants) are mapped and fitted; the
  AA references in `AA-reference/` point to the functional motor trajectories.
- The `rotate/` subfolder contains the CG/AA rotation experiments described in
  the original README; use `rotate_CG.sh` there for the validated CG rotation
  protocol.  AA rotation scripts remain experimental.
- The `gyration/` subfolder supplements the SASA measurements with AA vs CG
  radius-of-gyration comparisons (`analysis.ipynb`).
