"""
Reads file inp and writes to file out, changing all values in column "mass" (index 8) of the [ atoms ] section to 1.0.
"""

import sys
import re
import os
inp = sys.argv[1]
out = sys.argv[2]

CHECK=0
with open(inp, 'r') as file, open(out, 'w+') as fileout:
    for line in file:
        #Check if we have reached the [ atoms ] section
        if "[ atoms ]" in line:
            CHECK=1
            fileout.write(line)
        elif "[ bonds ]" in line:
            CHECK=0
            fileout.write(line)
        #If we are in the [ atoms ] section, if first column is a number, then change the mass to 1.0 and write the line
        elif CHECK == 1:
            data = re.split('\s+', line)
            if data[1].isdigit():
                data[8] = "1.0"
                fileout.write(" ".join(data) + "\n")
            else:
                fileout.write(line)
        else:
            fileout.write(line)