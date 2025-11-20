#!/bin/bash

# Define search directory
TRAJ_DIR="../trajectories"

# Find all directories in the trajectories folder that contain both a PDB and XTC file
for dir in "$TRAJ_DIR"/*/ ; do
    # Extract directory name without path
    dirname=$(basename "$dir")
    
    # Check if both required files exist
    if [ -f "$dir/$dirname.pdb" ] && [ -f "$dir/$dirname.xtc" ]; then
        echo "Processing $dirname..."
        
        # Run the Python script with proper arguments
        python measure_parameters.py "$dir/$dirname.pdb" "$dir/$dirname.xtc" "$dirname" $1
        
        # Check if command was successful
        if [ $? -eq 0 ]; then
            echo "✓ Successfully processed $dirname"
        else
            echo "✗ Error processing $dirname (exit code: $?)"
        fi
    else
        echo "⚠ Skipping $dirname - missing required files"
    fi
done

echo "All processing complete!"