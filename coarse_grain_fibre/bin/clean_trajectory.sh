#!/bin/bash

# Script to process GROMACS trajectory files
# - Uses index file for reference
# - Makes all molecules whole in all frames
# - Centers molecules in xy-plane based on COM of non-counterion molecules

# Exit on error
set -e

# Function to display usage
usage() {
    echo "Usage: $0 -s <structure.pdb> -f <trajectory.xtc> -p <topology.top> -t <topology.tpr> -o <output_name>"
    echo "Options:"
    echo "  -s: Structure file (.pdb)"
    echo "  -f: Trajectory file (.xtc)"
    echo "  -p: Topology file (.top)"
    echo "  -t: Run input file (.tpr)"
    echo "  -o: Output name (files will be saved in ../trajectories/OUTPUT_NAME/)"
    echo "  -h: Display this help message"
    exit 1
}

# Parse command line arguments
while getopts "s:f:p:t:n:o:h" opt; do
    case ${opt} in
        s )
            STRUCTURE=$OPTARG
            ;;
        f )
            TRAJ=$OPTARG
            ;;
        p )
            TOP=$OPTARG
            ;;
        t )
            TPR=$OPTARG
            ;;
        o )
            OUTPUT_NAME=$OPTARG
            ;;
        h )
            usage
            ;;
        \? )
            usage
            ;;
    esac
done

# Check if required parameters are provided
if [ -z "$STRUCTURE" ] || [ -z "$TRAJ" ] || [ -z "$TOP" ] || [ -z "$TPR" ] || [ -z "$OUTPUT_NAME" ]; then
    echo "Error: Missing required parameters" 
    echo $STRUCTURE $TRAJ $TOP $TPR $OUTPUT_NAME
    usage
fi

# Check if input files exist
for file in "$STRUCTURE" "$TRAJ" "$TOP" "$TPR"; do
    if [ ! -f "$file" ]; then
        echo "Error: File $file does not exist"
        exit 1
    fi
done

# Create output directory
OUTPUT_DIR="../trajectories/${OUTPUT_NAME}"
NATOMS=28
mkdir -p "$OUTPUT_DIR"
mkdir -p ../logs
rm -f ../tmp/*#*

echo "Processing trajectory..."
echo "Input structure: $STRUCTURE"
echo "Input trajectory: $TRAJ"
echo "Input topology: $TOP"
echo "Input TPR: $TPR"
echo "Output directory: $OUTPUT_DIR"
NDX=../CG-reference/index.ndx
# Step 1: Create a temporary directory for intermediate files 
TEMP_DIR=../tmp
printf "11\n" | gmx editconf -f $STRUCTURE -n $NDX -o ${TEMP_DIR}/${OUTPUT_NAME}.pdb
STRUCTURE=${TEMP_DIR}/${OUTPUT_NAME}.pdb
cp $TRAJ ${TEMP_DIR}/${OUTPUT_NAME}.xtc
TRAJ=${TEMP_DIR}/${OUTPUT_NAME}.xtc
cp $TOP ${TEMP_DIR}/${OUTPUT_NAME}.top
TOP=${TEMP_DIR}/${OUTPUT_NAME}.top
cp $TPR ${TEMP_DIR}/${OUTPUT_NAME}.tpr
TPR=${TEMP_DIR}/${OUTPUT_NAME}.tpr

gmx trjconv -f "$TRAJ" -s "$STRUCTURE" -pbc nojump -o ${OUTPUT_DIR}/${OUTPUT_NAME}_nojump.xtc  << EOF
0
EOF

# First remove PBC jumps
echo "Removing jumps across periodic boundaries..."
gmx trjconv -f "$TRAJ" -s "$STRUCTURE" -pbc whole -o ${TEMP_DIR}/whole.xtc << EOF
0
EOF

echo "Removing jumps across periodic boundaries..."
gmx trjconv -f ${TEMP_DIR}/whole.xtc -s "$TPR" -n $NDX -pbc cluster -center -o ${TEMP_DIR}/cluster.xtc << EOF
11
11
11
EOF

gmx trjconv -f ${TEMP_DIR}/cluster.xtc -s "$TPR" -pbc mol -n $NDX -o ${OUTPUT_DIR}/${OUTPUT_NAME}.xtc  << EOF
11
EOF

declare -A molecule_resid_map  # Maps molecule name to its assigned resid
current_mol_idx=1             # Counter for unique molecule types

# Build atom-to-resid mapping
declare -A atom_to_resid       # Maps atom index to its new resid
current_atom=1                # Tracks current atom position

# Find and process the molecules section in TOP file
found_molecules=false
while IFS= read -r line; do
    if ! $found_molecules; then
        if [[ $line == *"[ molecules ]"* ]]; then
            found_molecules=true
        fi
        continue
    fi
    
    # Exit if we hit another section
    if [[ $line == *"["* && $line != *"[ molecules ]"* ]]; then
        break
    fi
    
    # Skip empty lines and comments
    if [[ -z $line || $line == \;* ]]; then
        continue
    fi
    
    # Extract molecule name and count
    mol_name=""
    mol_count=""
    read -r mol_name mol_count <<< "$line"
    
    # Skip if not a valid molecule line
    if [[ -z $mol_name || -z $mol_count ]]; then
        continue
    fi
    
    # Assign resid to this molecule type if not already assigned
    if [[ -z ${molecule_resid_map[$mol_name]} ]]; then
        molecule_resid_map[$mol_name]=$current_mol_idx
        current_mol_idx=$((current_mol_idx + 1))
        echo "Molecule type: $mol_name -> ResID: ${molecule_resid_map[$mol_name]}"
    fi
    
    # Get the resid for this molecule type
    mol_resid=${molecule_resid_map[$mol_name]}
    
    # Calculate atoms in this molecule group
    atoms_in_group=0
    if [[ $mol_name == "NA" ]]; then
        atoms_in_group=$mol_count  # 1 atom per NA molecule
    else
        atoms_in_group=$((NATOMS * mol_count))
    fi
    
    # Assign the same resid to all atoms in this block
    for ((i=0; i<atoms_in_group; i++)); do
        atom_idx=$((current_atom + i))
        atom_to_resid[$atom_idx]=$mol_resid
    done
    
    # Update current atom counter
    current_atom=$((current_atom + atoms_in_group))
    
    echo "Processed: $mol_name ($mol_count), Atom range: $((current_atom-atoms_in_group))-$((current_atom-1)), ResID: $mol_resid"
    
done < "$TOP"

total_atoms=$((current_atom - 1))
echo "Total atoms mapped: $total_atoms"

# Process the PDB file and update residue numbers
atom_counter=1
TMP_PDB=../tmp/tmp.pdb
rm -f $TMP_PDB
while IFS= read -r line; do
    if [[ ${line:0:4} == "ATOM" ]]; then
        # Get the resid for this atom
        new_resid=${atom_to_resid[$atom_counter]}
        
        # If no mapping found, keep original
        if [[ -z $new_resid ]]; then
            echo "${line}"
            printf "${line}\n" >> "$TMP_PDB"
        else
            # Format the new residue number (right-justified in 4 spaces)
            formatted_resid=$(printf "%4d" $new_resid)
            
            # Replace the residue number in the line
            new_line="${line:0:22}${formatted_resid}${line:26}"
            printf "${new_line}\n" >> "$TMP_PDB"
        fi
        
        # Increment atom counter
        atom_counter=$((atom_counter + 1))
    else
        # Keep non-atom lines as they are
        printf "${line}\n" >> "$TMP_PDB"
    fi
done < $STRUCTURE

# Move the temporary file to the output file
cp "$TMP_PDB" ${OUTPUT_DIR}/${OUTPUT_NAME}.pdb

cp $NDX ${OUTPUT_DIR}/${OUTPUT_NAME}.ndx
cp $TOP ${OUTPUT_DIR}/${OUTPUT_NAME}.top
rm -f "$TEMP_INDICES"
rm -f ${TEMP_DIR}/${OUTPUT_NAME}_renamed.pdb


echo "Done!"
echo "Processed files are available in $OUTPUT_DIR"
echo "  - Structure: ${OUTPUT_DIR}/${OUTPUT_NAME}.pdb"
echo "  - Trajectory: ${OUTPUT_DIR}/${OUTPUT_NAME}.xtc"
echo "  - Topology: ${OUTPUT_DIR}/${OUTPUT_NAME}.top"
echo "  - Index: ${OUTPUT_DIR}/${OUTPUT_NAME}.ndx"