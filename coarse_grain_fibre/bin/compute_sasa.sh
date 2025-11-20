#!/usr/bin/env bash
set -euo pipefail
# NOT USED IN PRODUCTION BUT SHOULD BE USED
# Usage: ./run_sasa.sh STRUCT(PDB/GRO) TRAJ_XTC OUTPUT_BASENAME [DT_ps]
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 structure.(pdb|gro) trajectory.xtc output_basename [dt_ps]"
  exit 1
fi

STRUCT_IN="$1"
TRAJ_IN="$2"
OUTBASE="$3"
DT_PS="${4:-200000}"   # default: 200000 ps

# Tunables
PROBE=0.185
NDOTS=4800
NBOX_X=1
NBOX_Y=1
NBOX_Z=3

# Paths
SASA_DIR="../SASA"
mkdir -p "${SASA_DIR}"

RADII_DEST="./vdwradii.dat"
TRAJ_THIN="${SASA_DIR}/${OUTBASE}_thin.xtc"
NDX_FILE="${SASA_DIR}/${OUTBASE}_mid.ndx"
SASA_OUT="${SASA_DIR}/${OUTBASE}.xvg"

echo ">>> Copying vdW radii table"
if [[ -f ./vdwradii.dat ]]; then
  echo "./vdwradii.dat Exists"
elif [[ -f ../structure_files/vdwradii_CG.dat  ]]; then
  cp ../structure_files/vdwradii_CG.dat "${RADII_DEST}"
  echo "Copied ../structure_files/vdwradii_CG.dat  -> ${RADII_DEST}"
else
  echo "WARNING: No vdwradii.dat beside script or at ../vdwradii_CG.dat."
  echo "         Proceeding with GROMACS defaults."
fi

echo ">>> Thinning trajectory to every ${DT_PS} ps"
# group 0 = System
gmx trjconv -s "${STRUCT_IN}" -f "${TRAJ_IN}" -o "${TRAJ_THIN}" -dt "${DT_PS}" -tu ps <<< "0"

echo ">>> Collecting times from thinned trajectory"
# Extract frame times (ps) from the thinned traj
# Works across GROMACS versions by asking trjconv to write times and grepping them.
TMP_TIMES="${SASA_DIR}/${OUTBASE}_times.txt"
gmx trjconv -s "${STRUCT_IN}" -f "${TRAJ_THIN}" -o /dev/null -dump -1 2>/dev/null <<< "0" || true
# Use gmx check to print time stamps
gmx check -f "${TRAJ_THIN}" 2>/dev/null | awk '
  /t =/ { 
    # lines like: "  time  200000 (ps)  step ..."
    # fallback: catch "t = 2e+05"
  }
' >/dev/null || true

# Portable way: print one small file per frame to capture time with -b/-e
# We’ll iterate frames by index and query the time from gmx trjconv output each loop.

echo "@    title \"SASA vs time\"" > "${SASA_OUT}"
echo "@    xaxis  label \"Time (ps)\"" >> "${SASA_OUT}"
echo "@    yaxis  label \"SASA (nm^2)\"" >> "${SASA_OUT}"

FRAME_IDX=0
FIRST_DONE=0

# Ask gmx check for the number of frames
NFRAMES=$(gmx check -f "${TRAJ_THIN}" 2>/dev/null | awk '/Step +[0-9]+/ {last=$2} END{print (last==""?0:last+1)}')
if [[ "${NFRAMES}" -eq 0 ]]; then
  # Fallback method to count frames with trjconv -sep
  TMP_SEP_DIR="${SASA_DIR}/_sep_$$"
  mkdir -p "${TMP_SEP_DIR}"
  gmx trjconv -s "${STRUCT_IN}" -f "${TRAJ_THIN}" -o "${TMP_SEP_DIR}/frame.gro" -sep <<< "0"
  NFRAMES=$(ls -1 "${TMP_SEP_DIR}"/frame*.gro 2>/dev/null | wc -l | tr -d ' ')
  rm -rf "${TMP_SEP_DIR}"
fi

echo ">>> Found ${NFRAMES} frames in thinned trajectory"

while [[ "${FRAME_IDX}" -lt "${NFRAMES}" ]]; do
  echo ">>> Processing frame ${FRAME_IDX}/${NFRAMES}"

  # Extract single frame by index using -dump from gmx trjconv: get its time first
  # We get time by asking trjconv to dump that frame and parse stdout
  # Trick: use -skip/ -tu ps to position; but -dump expects time, not index.
  # So we derive time = FRAME_IDX * DT_PS relative to first frame.
  TIME_PS=$(( FRAME_IDX * DT_PS ))

  FRAME_GRO="${SASA_DIR}/_frame_${FRAME_IDX}.gro"
  REPL_GRO="${SASA_DIR}/_frame_${FRAME_IDX}_1x1x3.gro"
  TMP_XVG="${SASA_DIR}/_sasa_${FRAME_IDX}.xvg"

  # Extract exactly this frame from the thinned traj
  gmx trjconv -s "${STRUCT_IN}" -f "${TRAJ_THIN}" -o "${FRAME_GRO}" -b "${TIME_PS}" -e "${TIME_PS}" -tu ps <<< "0"

  # Replicate 1x1x3 on this single-frame structure
  gmx genconf -f "${FRAME_GRO}" -o "${REPL_GRO}" -nbox ${NBOX_X} ${NBOX_Y} ${NBOX_Z}

  # Build index once (first frame only)
  if [[ "${FIRST_DONE}" -eq 0 ]]; then
    echo ">>> Creating index (System; MID = resnr 5601..11200)"
    printf "a 5601-11200\nname 6 MID\nq\n" | gmx make_ndx -f "${REPL_GRO}" -o ${NDX_FILE}
    FIRST_DONE=1
  fi

  # Run SASA for this single frame; run in SASA_DIR so vdwradii.dat is picked up
  pushd "${SASA_DIR}" >/dev/null
  printf "0\n6\n" | gmx sasa \
    -s "$(basename "${REPL_GRO}")" \
    -f "$(basename "${REPL_GRO}")" \
    -n "$(basename "${NDX_FILE}")" \
    -ndots ${NDOTS} \
    -probe ${PROBE} \
    -o "$(basename "${TMP_XVG}")" \
    >/dev/null
  popd >/dev/null

  # Parse total SASA from the xvg (last non-comment line)
  SASA_VAL=$(awk 'BEGIN{val=""} /^[^@#]/ {val=$NF} END{print val}' "${TMP_XVG}")

  # Append time and SASA to final xvg
  printf "%d %.6f\n" "${TIME_PS}" "${SASA_VAL}" >> "${SASA_OUT}"

  # Cleanup temps for this frame
  rm -f "${FRAME_GRO}" "${REPL_GRO}" "${TMP_XVG}"

  FRAME_IDX=$(( FRAME_IDX + 1 ))
done

rm $RADII_DEST
echo ">>> Done. SASA timeseries written to: ${SASA_OUT}"
