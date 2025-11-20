#Reads input gro-file and topology file and reorders gro and topology file so that all groups in the topology file are continuous
#Usage: python reorder_gro_file.py -f <gro-file> -p <topology-file> -o <output-gro> -po <output-topology>

import sys, argparse, re

#Read inputs
parser = argparse.ArgumentParser()
parser.add_argument("-f","--gro_file", help="Input gro-file")
parser.add_argument("-p", "--topology", help="Input topology file")
parser.add_argument("-o", "--output_gro", help="Output gro-file")
parser.add_argument("-po", "--output_topology", help="Output topology file")
args = parser.parse_args()

MOLS = {'MOTOR-R':28, 'MOTOR-S':28,'NA':1,'W':1}

#Read gro-file
with open(args.gro_file,'r') as f, open(args.output_gro, 'w') as out_gro, open(args.output_topology, 'w') as out_top, open(args.topology, 'r') as top:
    #Read all lines in topology file after line including "[ molecules ]" and save as list
    top_lines = top.readlines()
    
    #Write all lines in topology file before line including "[ molecules ]" to output topology file
    out_top.writelines(top_lines[:top_lines.index("[ molecules ]\n")+2])
    top_lines = top_lines[top_lines.index("[ molecules ]\n")+2:]
    #Read all lines in gro-file and save as list
    gro_lines = f.readlines()

    #Write first two lines of gro-file to output gro-file
    out_gro.write(gro_lines[0])
    out_gro.write(gro_lines[1])

    #For each moleculetype in MOLS, find all atoms in gro-file and topology file and write to output gro-file and topology file
    OUTINDEX = 1
    RESNR = 1
    for mol,_ in MOLS.items():
        #Find position of atoms in gro-file from topology file
        INDEX = 2
        NMOLS = 0
        for line in top_lines:
            moltype, num = re.split(r'\s+',line.strip())
            num = int(num)
            if mol == moltype:
                #Find position of atoms in gro-file and write to output gro-file, reindexing atoms
                for i in range(INDEX, INDEX+num*MOLS[moltype]):
                    out_gro.write(str(RESNR).rjust(5)+gro_lines[i][5:15]+str(OUTINDEX).rjust(5)+gro_lines[i][20:])
                    OUTINDEX += 1
                    if OUTINDEX % 100000 == 0:
                        OUTINDEX = 0
                    if (i - INDEX + 1) % MOLS[moltype] == 0:
                        RESNR += 1
                NMOLS += num
            INDEX += num*MOLS[moltype]
        #Write number of molecules to output topology file
        out_top.write(mol.ljust(10)+str(NMOLS).rjust(5)+'\n')
    #Write last line of gro-file to output gro-file
    out_gro.write(gro_lines[-1])
                

