#!/bin/bash

set -euo pipefail

TRAJ_DIR="../trajectories"

for dir in "$TRAJ_DIR"/*/ ; do
    dirname=$(basename "$dir")

    pdb="$dir/$dirname.pdb"
    xtc="$dir/$dirname.xtc"

    if [[ -f "$pdb" && -f "$xtc" ]]; then
        echo "Processing $dirname..."
        python measure_parameters.py "$pdb" "$xtc" "$dirname" ALL
        echo "✓ Finished $dirname"
    else
        echo "⚠ Skipping $dirname - missing $dirname.pdb or $dirname.xtc"
    fi
done

echo "All processing complete!"
