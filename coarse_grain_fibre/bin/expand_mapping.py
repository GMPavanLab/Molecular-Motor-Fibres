# Purpose: Expand mapping file from one molecule to many molecules. Molecules are numbered from 1 to natoms and come two  blocks of 110 molecules, separated by a gap of 220 atoms.
import sys
import re

inp = sys.argv[1]
out = sys.argv[2]
natoms = int(sys.argv[3])
edit = 110
if len(sys.argv) > 4:
    edit = int(sys.argv[4])
print(edit)

with open(inp, 'r') as file, open(out, 'w+') as fileout:
    lines = file.readlines()
    for i in range(edit):
        for line in lines:
            #If line does not contain digits, then write it
            if (not any(char.isdigit() for char in line)) or '[' in line:
                if '[' in line:
                    name = line.strip('[]\n')
            else:
                #Extract all digits from line
                data = [int(x) for x in re.findall(r'\d+', line)]
                fileout.write(f'[{name}{str(i+1)}]\n')
                fileout.write(' '.join([str(x + natoms*i) for x in data]) + '\n')
                fileout.write('\n')
    for i in range(edit):
        for line in lines:
            #If line does not contain digits, then write it
            if not any(char.isdigit() for char in line) or '[' in line:
                if '[' in line:
                    name = line.strip('[]\n')
            else:
                data = [int(x) for x in re.findall(r'\d+', line)]
                fileout.write(f'[{name}{str(i+1+edit)}]\n')
                fileout.write(' '.join([str(x + edit*2 + natoms*(edit + i)) for x in data]) + '\n')
                fileout.write('\n')