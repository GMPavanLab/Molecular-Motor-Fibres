"""
Adds include statement to bottom of topol file. 
Usage:
python include_itp_in_top.py file.itp file.top n_molecules NAME_OF_RESTRAINT-file opt:MOLNAME"""
import sys
fileitp = sys.argv[1]
filetop = sys.argv[2]
n_molecules = int(sys.argv[4])
RESTRAINT_NAME = sys.argv[3]

resstring = f'\n\
#ifdef {RESTRAINT_NAME}\n\
#include "../structure_files/{RESTRAINT_NAME}.itp"\n\
#endif\n\n'

if len(sys.argv) == 6:
    mol_name = sys.argv[5]
else:
    with open(fileitp) as itp:
        mol_name = itp.readlines()[2].split()[0]

with open(filetop, 'r') as top:
    lines = top.readlines()

with open(filetop, 'w+') as top:
    system_idx = lines.index('[ system ]\n')
    if RESTRAINT_NAME != 'no':
        lines.insert(system_idx, resstring)
    lines.insert(system_idx, f'#include "{fileitp}"\n')
    lines.append(f"{mol_name}   {n_molecules}\n")
    top.writelines(lines)
