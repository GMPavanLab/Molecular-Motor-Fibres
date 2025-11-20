import sys
import re

filename = sys.argv[1]
residue_nmbr = 0
residue_map = {}
atom_residue_map = {}
check = False
#Iterate over lines in file with filename
with open(filename, 'r') as file:
    for line in file:
        if "@<TRIPOS>BOND" in line:
            check=False
        
        if check and len(re.split('\s+', line)) > 5:
            params = re.split('\s+', line)
            element = params[2][0]
            residue = params[8]
            if not residue in residue_map:
                residue_map[residue] = residue_nmbr
                atom_residue_map[residue] = {}
                residue_nmbr += 100
            if not element in atom_residue_map[residue]:
                atom_residue_map[residue][element] = 1
            if not element == "H":
                params[2] = params[2][0] + str(residue_map[residue] + atom_residue_map[residue][element])
            else:
                params[2] = params[2][0]
            atom_residue_map[residue][element] += 1
            
            print('  '.join(params))

        else:
            print(line, end='')

        if "@<TRIPOS>ATOM" in line:
            check = True
