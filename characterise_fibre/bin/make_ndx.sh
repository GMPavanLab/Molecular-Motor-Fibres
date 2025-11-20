#!/bin/bash

# Check if input files are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <topology_file> <pdb_file> <output_ndx_file> [natoms]"
    echo "  topology_file: GROMACS topology file (.top)"
    echo "  pdb_file: GROMACS PDB file (.pdb)"
    echo "  output_ndx_file: Output index file (.ndx)"
    echo "  natoms: Number of atoms in each motor molecule (optional)"
    exit 1
fi
TOP_FILE="$1"
PDB_FILE="$2"
OUTPUT_NDX="$3"
NATOMS=${4:-152}  # Default value if not provided

rm ../tmp/\#* 2>/dev/null
rm ../tmp/tmp* 2>/dev/null

# Initialize atom index counter and string variables for each group
ATOM_INDEX=1
MOTOR=""
MOTOR2M=""
MOTORB=""
ISOR=""
ISOS=""
LLG=""
RLG=""
HED=""

# Extract molecule information from topology file
echo "Parsing topology file: $TOP_FILE"
# Skip to the [ molecules ] section
MOLECULES=$(awk '/\[ molecules \]/,/^$/' "$TOP_FILE" | grep -v '^\s*;' | grep -v '^\s*$' | grep -v '\[ molecules \]')

# Initialize variables to store unique compound bases (without -R/-S)
declare -A COMPOUND_BASES
declare -A COMPOUND_COUNTS

# First pass: Identify unique compound bases
while read -r line; do
    # Skip comments and empty lines
    [[ "$line" =~ ^[[:space:]]*\; ]] && continue
    [[ -z "$line" ]] && continue
    
    COMPOUND=$(echo "$line" | awk '{print $1}')
    COUNT=$(echo "$line" | awk '{print $2}')
    
    if [[ "$COMPOUND" == *"motor"* ]]; then
        # Extract base compound name (without -R/-S)
        BASE_COMPOUND=$(echo "$COMPOUND" | sed 's/-[RS]$//')
        COMPOUND_BASES["$BASE_COMPOUND"]=1
        
        # Track total count for each compound
        if [[ -z "${COMPOUND_COUNTS[$BASE_COMPOUND]}" ]]; then
            COMPOUND_COUNTS["$BASE_COMPOUND"]=$COUNT
        else
            COMPOUND_COUNTS["$BASE_COMPOUND"]=$((${COMPOUND_COUNTS["$BASE_COMPOUND"]} + $COUNT))
        fi
    fi
done <<< "$MOLECULES"

# Get the unique compound bases
UNIQUE_BASES=()
for base in "${!COMPOUND_BASES[@]}"; do
    UNIQUE_BASES+=("$base")
done

echo "Found unique motor compound bases: ${UNIQUE_BASES[*]}"

# Process molecules and build atom selection strings
while read -r line; do
    # Skip comments and empty lines
    [[ "$line" =~ ^[[:space:]]*\; ]] && continue
    [[ -z "$line" ]] && continue
    
    COMPOUND=$(echo "$line" | awk '{print $1}')
    COUNT=$(echo "$line" | awk '{print $2}')
    
    # Process only motor-containing compounds
    if [[ "$COMPOUND" == *"motor"* ]]; then
        echo "Processing compound: $COMPOUND with $COUNT molecules"
        
        # Extract base compound name (without -R/-S)
        BASE_COMPOUND=$(echo "$COMPOUND" | sed 's/-[RS]$//')
        
        # Process each molecule of this compound
        for ((i=0; i<$COUNT; i++)); do
            START_ATOM=$ATOM_INDEX
            END_ATOM=$((ATOM_INDEX + NATOMS - 1))
            
            # Add atoms to MOTOR group
            if [[ -n "$MOTOR" ]]; then
                MOTOR="$MOTOR | a $START_ATOM-$END_ATOM"
            else
                MOTOR="a $START_ATOM-$END_ATOM"
            fi
            
            # Add atoms to the specific base compound group
            for base in "${UNIQUE_BASES[@]}"; do
                if [[ "$BASE_COMPOUND" == "$base" ]]; then
                    if [[ -n "${!base}" ]]; then
                        eval "$base=\"\$$base | a $START_ATOM-$END_ATOM\""
                    else
                        eval "$base=\"a $START_ATOM-$END_ATOM\""
                    fi
                    break
                fi
            done
            
            # Add atoms to isomer groups
            if [[ "$COMPOUND" == *"-R" ]]; then
                if [[ -n "$ISOR" ]]; then
                    ISOR="$ISOR | a $START_ATOM-$END_ATOM"
                else
                    ISOR="a $START_ATOM-$END_ATOM"
                fi
            elif [[ "$COMPOUND" == *"-S" ]]; then
                if [[ -n "$ISOS" ]]; then
                    ISOS="$ISOS | a $START_ATOM-$END_ATOM"
                else
                    ISOS="a $START_ATOM-$END_ATOM"
                fi
            fi
            
            # Add atoms to LLG/RLG/HED groups using per-molecule ranges (same style)
            # LLG: indices 72-105 within each molecule => +73 .. +106 (1-based inclusive)
            LLG_ATOMS="a $((START_ATOM + 72))-$((START_ATOM + 105))"
            # RLG: indices 118-151 within each molecule => +119 .. +152
            RLG_ATOMS="a $((START_ATOM + 118))-$((START_ATOM + 151))"
            # HED: indices 0-37 within each molecule => +1 .. +38
            HED_ATOMS="a $((START_ATOM + 0))-$((START_ATOM + 37))"

            if [[ -n "$LLG" ]]; then
                LLG="$LLG | $LLG_ATOMS"
            else
                LLG="$LLG_ATOMS"
            fi

            if [[ -n "$RLG" ]]; then
                RLG="$RLG | $RLG_ATOMS"
            else
                RLG="$RLG_ATOMS"
            fi

            if [[ -n "$HED" ]]; then
                HED="$HED | $HED_ATOMS"
            else
                HED="$HED_ATOMS"
            fi
            
            # Update atom index for next molecule
            ATOM_INDEX=$((END_ATOM + 1))
        done
    else
        # Skip non-motor compounds but update atom index
        # SOL typically has 3 atoms, NA has 1
        if [[ "$COMPOUND" == "SOL" ]]; then
            ATOM_INDEX=$((ATOM_INDEX + COUNT * 3))
        elif [[ "$COMPOUND" == "NA" ]]; then
            ATOM_INDEX=$((ATOM_INDEX + COUNT))
        else
            # Default assumption: 1 atom per molecule
            ATOM_INDEX=$((ATOM_INDEX + COUNT))
            echo "Warning: Unknown compound $COMPOUND, assuming 1 atom per molecule"
        fi
    fi
done <<< "$MOLECULES"

# Create temporary input file for gmx make_ndx
TMP_INPUT=../tmp/tmp1.txt
NDX=../tmp/tmp_ndx.ndx

# First make a copy of the existing groups
printf "del 7-10\n" > "$TMP_INPUT"
printf "0 | 1 | 2 | 3 | 4 | 5 | 6\n" >> "$TMP_INPUT"
printf "name 7 ALL\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp2.txt
# Add MOTOR group
printf "$MOTOR \n" >> "$TMP_INPUT"
printf "name 8 MOTOR\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

# Add specific compound groups
group_num=9
for base in "${UNIQUE_BASES[@]}"; do
    TMP_INPUT=../tmp/tmp_${base}.txt
    printf "${!base} \n" >> "$TMP_INPUT"
    printf "name $group_num $base \nq\n" >> "$TMP_INPUT"
    group_num=$((group_num + 1))
    gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
    cp "$OUTPUT_NDX" "$NDX"
done

TMP_INPUT=../tmp/tmp_4.txt

# Add isomer groups
printf "$ISOR \n" >> "$TMP_INPUT"
printf "name $group_num iso-R\nq\n" >> "$TMP_INPUT"
group_num=$((group_num + 1))
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_5.txt

printf "$ISOS \n" >> "$TMP_INPUT"
printf "name $group_num iso-S\nq\n" >> "$TMP_INPUT"
group_num=$((group_num + 1))
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

###############################################################################
# LLG: split string into three parts and build final LLG (unchanged pattern)
###############################################################################
# Split LLG string into three parts
LLG_length=${#LLG}
split_point1=$((LLG_length / 3))
split_point2=$((2 * LLG_length / 3))

clean_split1=$(echo "$LLG" | grep -b -o "| a" | awk -F: -v sp=$split_point1 '{if ($1 > sp) {print $1; exit}}')
if [ -n "$clean_split1" ]; then
    LLG_part1="${LLG:0:$clean_split1}"
    clean_split1=$((clean_split1 + 2))
else
    LLG_part1="${LLG:0:$split_point1}"
    clean_split1=$split_point1
fi

clean_split2=$(echo "$LLG" | grep -b -o "| a" | awk -F: -v sp=$split_point2 '{if ($1 > sp) {print $1; exit}}')
if [ -n "$clean_split2" ]; then
    LLG_part2="${LLG:$clean_split1:$((clean_split2-clean_split1))}"
    clean_split2=$((clean_split2 + 2))
else
    LLG_part2="${LLG:$clean_split1:$((split_point2-clean_split1))}"
    clean_split2=$split_point2
fi
LLG_part3="${LLG:$clean_split2}"

# Process first part of LLG
TMP_INPUT=../tmp/tmp_81.txt
printf "$LLG_part1 \n" >> "$TMP_INPUT"
printf "name $group_num LLG_part1\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

# Process second part of LLG
TMP_INPUT=../tmp/tmp_82.txt
printf "$LLG_part2 \n" > "$TMP_INPUT"
printf "name $((group_num + 1)) LLG_part2\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

# Process third part of LLG
TMP_INPUT=../tmp/tmp_83.txt
printf "$LLG_part3 \n" > "$TMP_INPUT"
printf "name $((group_num + 2)) LLG_part3\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

# Join the three parts to create the final LLG group
TMP_INPUT=../tmp/tmp_84.txt
printf "$((group_num)) | $((group_num + 1)) | $((group_num + 2))\n" > "$TMP_INPUT"
printf "name $((group_num + 3)) LLG\n" >> "$TMP_INPUT"
printf "del $((group_num))-$((group_num + 2))\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

# After creating final LLG, bump group_num (keep pattern)
group_num=$((group_num + 1))

###############################################################################
# RLG: follow LLG pattern exactly
###############################################################################
RLG_length=${#RLG}
split_point1=$((RLG_length / 3))
split_point2=$((2 * RLG_length / 3))

clean_split1=$(echo "$RLG" | grep -b -o "| a" | awk -F: -v sp=$split_point1 '{if ($1 > sp) {print $1; exit}}')
if [ -n "$clean_split1" ]; then
    RLG_part1="${RLG:0:$clean_split1}"
    clean_split1=$((clean_split1 + 2))
else
    RLG_part1="${RLG:0:$split_point1}"
    clean_split1=$split_point1
fi

clean_split2=$(echo "$RLG" | grep -b -o "| a" | awk -F: -v sp=$split_point2 '{if ($1 > sp) {print $1; exit}}')
if [ -n "$clean_split2" ]; then
    RLG_part2="${RLG:$clean_split1:$((clean_split2-clean_split1))}"
    clean_split2=$((clean_split2 + 2))
else
    RLG_part2="${RLG:$clean_split1:$((split_point2-clean_split1))}"
    clean_split2=$split_point2
fi
RLG_part3="${RLG:$clean_split2}"

TMP_INPUT=../tmp/tmp_91.txt
printf "$RLG_part1 \n" > "$TMP_INPUT"
printf "name $group_num RLG_part1\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_92.txt
printf "$RLG_part2 \n" > "$TMP_INPUT"
printf "name $((group_num + 1)) RLG_part2\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_93.txt
printf "$RLG_part3 \n" > "$TMP_INPUT"
printf "name $((group_num + 2)) RLG_part3\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_94.txt
printf "$((group_num)) | $((group_num + 1)) | $((group_num + 2))\n" > "$TMP_INPUT"
printf "name $((group_num + 3)) RLG\n" >> "$TMP_INPUT"
printf "del $((group_num))-$((group_num + 2))\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

group_num=$((group_num + 1))

###############################################################################
# HED: follow LLG pattern exactly
###############################################################################
HED_length=${#HED}
split_point1=$((HED_length / 3))
split_point2=$((2 * HED_length / 3))

clean_split1=$(echo "$HED" | grep -b -o "| a" | awk -F: -v sp=$split_point1 '{if ($1 > sp) {print $1; exit}}')
if [ -n "$clean_split1" ]; then
    HED_part1="${HED:0:$clean_split1}"
    clean_split1=$((clean_split1 + 2))
else
    HED_part1="${HED:0:$split_point1}"
    clean_split1=$split_point1
fi

clean_split2=$(echo "$HED" | grep -b -o "| a" | awk -F: -v sp=$split_point2 '{if ($1 > sp) {print $1; exit}}')
if [ -n "$clean_split2" ]; then
    HED_part2="${HED:$clean_split1:$((clean_split2-clean_split1))}"
    clean_split2=$((clean_split2 + 2))
else
    HED_part2="${HED:$clean_split1:$((split_point2-clean_split1))}"
    clean_split2=$split_point2
fi
HED_part3="${HED:$clean_split2}"

TMP_INPUT=../tmp/tmp_101.txt
printf "$HED_part1 \n" > "$TMP_INPUT"
printf "name $group_num HED_part1\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_102.txt
printf "$HED_part2 \n" > "$TMP_INPUT"
printf "name $((group_num + 1)) HED_part2\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_103.txt
printf "$HED_part3 \n" > "$TMP_INPUT"
printf "name $((group_num + 2)) HED_part3\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_104.txt
printf "$((group_num)) | $((group_num + 1)) | $((group_num + 2))\n" > "$TMP_INPUT"
printf "name $((group_num + 3)) HED\n" >> "$TMP_INPUT"
printf "del $((group_num))-$((group_num + 2))\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

group_num=$((group_num + 1))

###############################################################################
# COR = MOTOR \ (LLG ∪ RLG ∪ HED)
# MOTOR is group 8. The last three created groups (in order) are LLG, RLG, HED
# with indices (group_num-3), (group_num-2), (group_num-1) respectively.
###############################################################################
TMP_INPUT=../tmp/tmp_9.txt
printf "8 & ! $((group_num - 3)) & ! $((group_num - 2)) & ! $((group_num - 1))\n" >> "$TMP_INPUT"
printf "name $group_num COR \nq\n" >> "$TMP_INPUT"
group_num=$((group_num + 1))
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_10.txt

# Merge NA groups
printf "5\n" >> "$TMP_INPUT"
printf "name $group_num NA \nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"
cp "$OUTPUT_NDX" "$NDX"

TMP_INPUT=../tmp/tmp_11.txt
# Keep only the groups we want (delete the others)
printf "del 2-7\nq\n" >> "$TMP_INPUT"
gmx make_ndx -f "$PDB_FILE" -o "$OUTPUT_NDX" -n "$NDX" < "$TMP_INPUT"

# Run gmx make_ndx
echo "Running gmx make_ndx to generate the index file..."
echo "Index file successfully created: $OUTPUT_NDX"
