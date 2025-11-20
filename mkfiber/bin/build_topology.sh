#!/bin/bash
#UTILITY SCRIPT
#Builds topology file in the order of supplied topologies

TOPOLOGY_PATHS=()
NUMBER_MOLECULES=()
NAMES=()
NUMBER_IONS=0

while [[ $# -gt 0 ]]; do
  case $1 in
    -t|--topologies)
      shift #past argument
      while [[ $1 == *.top ]]; do
         TOPOLOGY_PATHS+=("$1")
         shift #past value
      done
      ;;
    -n|--number_of_molecules)
      shift #past argument
      while [[ $1 =~ ^[0-9]+$ ]]; do
         NUMBER_MOLECULES+=("$1")
         shift #past value
      done
      ;;
    -o |--output)
        shift
        OUTPUT_FILE_NAME=$1
        shift
        ;;
    -ion)
        shift
        NAME_ION=$1
        shift
        ;;
    -nc|--number_ions)
        shift
        NUMBER_IONS=$1
        shift
        ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
    *)
      echo "Unknown value $1"
      return 1
      ;;
  esac
done

for top in ${TOPOLOGY_PATHS[@]}; do
    x=${top%.top}
    NAMES+=("${x##*/}")
    OUTPATH=../structure_files/${x##*/}.itp
    python3 top2itp.py $top $OUTPATH
done

cat << EOF > $OUTPUT_FILE_NAME
;   This is a standalone topology file
;   Created by:   -, mkfiber/bin/build_topology.sh
;     -
;

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5          0.83333333  

[ atomtypes ]
; name    at.num    mass    charge ptype  sigma      epsilon
c3             6  12.010000  0.00000000  A     0.33977095      0.4510352
hc             1   1.008000  0.00000000  A      0.2600177      0.0870272
h1             1   1.008000  0.00000000  A     0.24219973      0.0870272
os             8  16.000000  0.00000000  A     0.31560978      0.3037584
ca             6  12.010000  0.00000000  A     0.33152123      0.4133792
ha             1   1.008000  0.00000000  A     0.26254785      0.0673624
ss            16  32.060000  0.00000000  A     0.35324134      1.1815616
ce             6  12.010000  0.00000000  A     0.33152123      0.4133792
cd             6  12.010000  0.00000000  A     0.33152123      0.4133792
c              6  12.010000  0.00000000  A     0.33152123      0.4133792
oh             8  16.000000  0.00000000  A     0.32428713       0.389112
ho             1   1.008000  0.00000000  A    0.053792465      0.0196648
o              8  16.000000  0.00000000  A     0.30481209      0.6121192
OW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01
HW           1       1.008   0.0000  A   0.00000e+00  0.00000e+00
Cl          17      35.45    0.0000  A   4.47766e-01  1.48913e-01
Na          11      22.99    0.0000  A   2.43928e-01  3.65846e-02
MG          12      24.305   0.0000  A   1.41225e-01  3.74342e+00
C0          20      40.08    0.0000  A   3.05240e-01  1.92376e+00
K           19      39.10    0.0000  A   3.03797e-01  8.10369e-01
Li           3       6.94    0.0000  A   2.02590e-01  7.65672e-02
Zn          30      65.4     0.0000  A   1.95998e-01  5.23000e-02


[ system ]
; Name
Motor fiber

[ molecules ]
; Compound       #mols
EOF

for i in "${!NAMES[@]}"; do
   gmx make_ndx -f ../../structure_files/pdb/${NAMES[$i]}.pdb -o ../tmp/index.ndx <<EOF
3 | 4
q
EOF
   echo 5 | gmx genrestr -f ../../structure_files/pdb/${NAMES[$i]}.pdb -o ../structure_files/POSRES_FIBRE.itp -fc 200 200 0 -n ../tmp/index.ndx
   python3 include_itp_in_top.py ../structure_files/${NAMES[$i]}.itp $OUTPUT_FILE_NAME POSRES_FIBRE ${NUMBER_MOLECULES[$i]}
done

if [ $NUMBER_IONS -gt 0 ]; then
    python3 include_itp_in_top.py ../structure_files/sol_ions_gaff.itp $OUTPUT_FILE_NAME no $NUMBER_IONS $NAME_ION
fi
