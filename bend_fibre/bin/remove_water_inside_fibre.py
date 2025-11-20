import sys
import re
import os
import numpy as np
inp = sys.argv[1]
out = sys.argv[2]
top = sys.argv[3]

NREM = 0
NATOMS=2000
XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN = -1, -1,-1, 10**9, 10**9,10**9

with open(inp, 'r') as file:
    for i, line in enumerate(file):
        if i == 1:
            NATOMS = int(line)
            print(f"{NATOMS} atoms to start")
        if i == NATOMS+2:
            X,Y,Z = [float(x) for x in re.split('\s+', line)[1:4]]
            print(X,Y,Z)
        if "MM2" in line:
            data = re.split('\s+', line)
            x,y,z = [float(x) for x in data[-4:-1]]
            #Update MIN and MAX
            XMAX = max(XMAX, x)
            YMAX = max(YMAX, y)
            ZMAX = max(ZMAX, z)
            XMIN = min(XMIN, x)
            YMIN = min(YMIN, y)
            ZMIN = min(ZMIN, z)

CHECK=0
margin=0.2
ZMAX = ZMAX -margin
ZMIN = ZMIN +margin
XMAX = XMAX -margin
XMIN = XMIN +margin
YMAX = YMAX -margin
YMIN = YMIN +margin

with open(inp, 'r') as file, open(out, 'w+') as fileout:
    for line in file:
        if "W" in line:
            data = re.split('\s+', line)
            x,y,z = [float(x) for x in data[-4:-1]]
            if ZMIN < z < ZMAX and XMIN < x < XMAX and YMIN < y < YMAX:
                NREM += 1
            else:
                fileout.write(line)
        else:
            fileout.write(line)

with open(out, "r") as file:
    with open(f"{out[:-3]}.txt", "w") as tmp_file:
        for i, line in enumerate(file):
            if i == 1:
                tmp_file.write(f"{NATOMS - NREM}\n")
            else:
                tmp_file.write(line)

os.replace(f"{out[:-3]}.txt", out)
print(f'{NREM} W removed')
CHECK=True

with open(top, "r") as file, open(f'{top[:-4]}tmp.top', "w") as tmp_file:
    for line in file:
        if CHECK and "W" in line:
            data = re.split('\s+', line)
            data[1] = str(int(data[1]) - NREM)
            print(f"{data[1]} W left")
            tmp_file.write(" ".join(data))
            CHECK=False
        else:
            tmp_file.write(line)

os.replace(top, f'{top[:-4]}bck.top')
os.replace(f'{top[:-4]}tmp.top', top)