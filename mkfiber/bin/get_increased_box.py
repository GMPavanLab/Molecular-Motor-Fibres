import sys
import re
path = sys.argv[1]

with open(path, 'r') as file:
    for line in file:
        if "CRYST" in line:
            clinc = [float(x) for x in re.split('\s+', line)[1:4]]
            clinc[0] = 12
            clinc[1] = 12
            clinc[2] *= 1/10
            print(' '.join([str(x) for x in clinc]))
            break