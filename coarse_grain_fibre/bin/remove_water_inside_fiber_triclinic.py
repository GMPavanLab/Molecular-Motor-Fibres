import sys
import re
import os
import numpy as np

inp = sys.argv[1]
out = sys.argv[2]
top = sys.argv[3]
r = float(sys.argv[4])
NREM = 0
NATOMS = 2000

# Function to calculate distance in triclinic coordinates (xy-plane only)
def distance_in_xy_plane(coords, box_center, box_matrix):
    # Get the coordinates in the xy-plane
    x, y = coords[0], coords[1]
    center_x, center_y = box_center[0], box_center[1]
    
    # Calculate the difference vector
    dx = x - center_x
    dy = y - center_y
    
    # Calculate fractional coordinates for the difference
    inv_box = np.linalg.inv(box_matrix)
    frac_diff = np.dot(inv_box, [dx, dy, 0])
    
    # Apply minimum image convention in fractional space
    frac_diff[0] = frac_diff[0] - np.round(frac_diff[0])
    frac_diff[1] = frac_diff[1] - np.round(frac_diff[1])
    
    # Convert back to Cartesian
    real_diff = np.dot(box_matrix, frac_diff)
    
    # Return the distance in the xy-plane
    return np.sqrt(real_diff[0]**2 + real_diff[1]**2)

# Read the file to get box dimensions and atom count
box_vectors = None
with open(inp, 'r') as file:
    for i, line in enumerate(file):
        if i == 1:
            NATOMS = int(line)
            print(f"{NATOMS} atoms to start")
        elif i == NATOMS + 2:
            # Get all box vectors for a triclinic box
            box_data = [float(x) for x in re.split(r'\s+', line.strip()) if x]
            
            # Properly parse the box vectors 
            # In GROMACS format, triclinic box is represented as:
            # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            # But sometimes only the first 9 or 6 values are present
            
            if len(box_data) >= 9:
                # Full triclinic box
                v1x, v2y, v3z = box_data[0], box_data[1], box_data[2]
                v1y, v1z = box_data[3], box_data[4]
                v2x, v2z = box_data[5], box_data[6]
                v3x, v3y = box_data[7], box_data[8]
            elif len(box_data) >= 6:
                # Partially specified triclinic box
                v1x, v2y, v3z = box_data[0], box_data[1], box_data[2]
                v1y, v1z = box_data[3], box_data[4]
                v2x = box_data[5] if len(box_data) > 5 else 0
                v2z = box_data[6] if len(box_data) > 6 else 0
                v3x = box_data[7] if len(box_data) > 7 else 0
                v3y = box_data[8] if len(box_data) > 8 else 0
            else:
                # Simple box, treat as rectangular
                v1x, v2y, v3z = box_data[0], box_data[1], box_data[2]
                v1y = v1z = v2x = v2z = v3x = v3y = 0
            
            # Create the box matrix
            box_matrix = np.array([
                [v1x, v2x, v3x],
                [v1y, v2y, v3y],
                [v1z, v2z, v3z]
            ])
            
            print(f"Box matrix:\n{box_matrix}")
            break

if box_matrix is None:
    print("Error: Could not extract box information correctly")
    sys.exit(1)

# The true center of a triclinic box is at (0.5, 0.5, 0.5) in fractional coordinates
# Convert this to Cartesian coordinates using the box matrix
frac_center = np.array([0.5, 0.5, 0.5])
box_center = np.dot(box_matrix, frac_center)
print(f"Box center: {box_center}")

# Process the file to remove waters in the cylinder
with open(inp, 'r') as file, open(out, 'w+') as fileout:
    for line in file:
        if "W" in line:
            data = re.split(r'\s+', line)
            if len(data) >= 4:  # Ensure there are enough elements
                coords = [float(x) for x in data[-4:-1]]
                
                # Calculate distance from center in the xy-plane only
                if distance_in_xy_plane(coords, box_center, box_matrix) < r:
                    NREM += 1
                else:
                    fileout.write(line)
            else:
                fileout.write(line)  # Write the line if it doesn't have proper coordinates
        else:
            fileout.write(line)

# Update atom count in output file
with open(out, "r") as file:
    with open(f"{out[:-3]}.txt", "w") as tmp_file:
        for i, line in enumerate(file):
            if i == 1:
                tmp_file.write(f"{NATOMS - NREM}\n")
            else:
                tmp_file.write(line)
os.replace(f"{out[:-3]}.txt", out)
print(f'{NREM} W removed')

# Update topology file
CHECK = True
with open(top, "r") as file, open(f'{top[:-4]}tmp.top', "w") as tmp_file:
    for line in file:
        if CHECK and "W" in line:
            data = re.split(r'\s+', line)
            data[1] = str(int(data[1]) - NREM)
            print(f"{data[1]} W left")
            tmp_file.write(" ".join(data))
            CHECK = False
        else:
            tmp_file.write(line)
os.replace(top, f'{top[:-4]}bck.top')
os.replace(f'{top[:-4]}tmp.top', top)