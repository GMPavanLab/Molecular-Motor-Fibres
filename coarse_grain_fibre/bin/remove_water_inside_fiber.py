import sys
import re
import os
inp = sys.argv[1]
out = sys.argv[2]
top = sys.argv[3]
r = float(sys.argv[4])

NREM = 0
NATOMS=2000

with open(inp, 'r') as file:
    for i, line in enumerate(file):
        if i == 1:
            NATOMS = int(line)
            print(f"{NATOMS} atoms to start")
        if i == NATOMS+2:
            X,Y,Z = [float(x) for x in re.split('\s+', line)[-4:-1]]
            print(X,Y,Z)

with open(inp, 'r') as file, open(out, 'w+') as fileout:
    for line in file:
        if "W" in line:
            data = re.split('\s+', line)
            x,y,z = [float(x) for x in data[-4:-1]]
            if ((x-X/2)**2 + (y-Y/2)**2)**0.5 < r:
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