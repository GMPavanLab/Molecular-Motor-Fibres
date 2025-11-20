import sys
filein = sys.argv[1]
fileout = sys.argv[2]
check = False
with open(filein, 'r') as inf, open(fileout, 'w') as outf:
    for line in inf:
        if '[ moleculetype ]' in line:
            check = True
        elif '[ system ]' in line:
            check = False
        if check:
            outf.write(line)
