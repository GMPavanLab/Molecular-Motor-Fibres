#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Check if arguments were provided
if [ $# -eq 0 ]; then
    echo "Error: No input FILENAMEs provided"
    echo "Usage: $0 FILENAME1 [FILENAME2 FILENAME3 ...]"
    exit 1
fi

# Process each FILENAME provided as an argument
for FILENAME in "$@"; do
    echo "Processing $FILENAME..."
    
    STRUCTURE="../../mkfiber/system/${FILENAME}/npt_${FILENAME}.pdb"
    TRAJ="../../mkfiber/system/${FILENAME}/npt_${FILENAME}.xtc"
    TOP="../../mkfiber/system/${FILENAME}/${FILENAME}.top"
    TPR="../../mkfiber/system/${FILENAME}/npt_${FILENAME}.tpr"
    
    # Check if required files exist
    if [ ! -f "$STRUCTURE" ]; then
        echo "Error: Structure file $STRUCTURE not found"
        exit 1
    fi
    
    if [ ! -f "$TRAJ" ]; then
        echo "Error: Trajectory file $TRAJ not found"
        exit 1
    fi
    
    if [ ! -f "$TOP" ]; then
        echo "Error: Topology file $TOP not found"
        exit 1
    fi
    
    if [ ! -f "$TPR" ]; then
        echo "Error: TPR file $TPR not found"
        exit 1
    fi
    
    # Run clean_trajectory.sh
    echo "Cleaning trajectory for $FILENAME..."
    bash clean_trajectory.sh -s "$STRUCTURE" -f "$TRAJ" -p "$TOP" -t "$TPR" -o "$FILENAME"
    
    # Set output directory (assuming OUTPUT_FILENAME should be FILENAME)
    OUTPUT_DIR="../trajectories/${FILENAME}"
    echo ${OUTPUT_DIR}/${FILENAME}.pdb
    
    # Run measure_parameters.py
    echo "Measuring parameters for $FILENAME..."
    python measure_parameters.py "${OUTPUT_DIR}/${FILENAME}.pdb" "${OUTPUT_DIR}/${FILENAME}.xtc" "$FILENAME" ALL
    
    echo "Completed processing $FILENAME"
    echo "----------------------------------------"
done

echo "All processing completed successfully!"