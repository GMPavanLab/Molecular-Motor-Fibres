#!/usr/bin/env python3
"""
Remove water molecules from within cylindrical fibers in GROMACS systems using MDAnalysis.
Output in PDB format.
"""

import numpy as np
import argparse
import sys
import re
import MDAnalysis as mda


def read_top_file(filename):
    """Read a GROMACS .top file and return its content."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    return lines


def update_water_count_in_top(lines, new_water_count):
    """Update the water molecule count in the topology file."""
    updated_lines = []
    in_molecules_section = False
    
    for line in lines:
        # Check if we're in the [ molecules ] section
        if '[ molecules ]' in line:
            in_molecules_section = True
            updated_lines.append(line)
            continue
        
        # If we're in the molecules section and find water (W)
        if in_molecules_section and not line.strip().startswith(';'):
            parts = line.split()
            if parts and parts[0] == 'W':
                # Update water count
                updated_lines.append(f"W    {new_water_count}\n")
                continue
            elif parts and parts[0].startswith('['):
                # We've left the molecules section
                in_molecules_section = False
        
        updated_lines.append(line)
    
    return updated_lines


def write_top_file(filename, lines):
    """Write the topology file."""
    with open(filename, 'w') as f:
        f.writelines(lines)


def calculate_cylinder_properties(atoms_group):
    """Calculate center of mass and z-extent for a cylinder."""
    # Calculate center of mass for x,y coordinates
    positions = atoms_group.positions
    com_x = np.mean(positions[:, 0])
    com_y = np.mean(positions[:, 1])
    
    # Find z-extent
    z_coords = positions[:, 2]
    z_min = np.min(z_coords)
    z_max = np.max(z_coords)
    
    return com_x, com_y, z_min, z_max


def is_water_inside_cylinder(water_pos, cylinder_com_x, cylinder_com_y, 
                           z_min, z_max, radius):
    """Check if a water atom is inside a cylinder."""
    # Check if z-coordinate is within cylinder height
    if water_pos[2] < z_min or water_pos[2] > z_max:
        return False
    
    # Calculate perpendicular distance to cylinder axis
    dx = water_pos[0] - cylinder_com_x
    dy = water_pos[1] - cylinder_com_y
    distance = np.sqrt(dx**2 + dy**2)
    
    return distance < radius


def main():
    parser = argparse.ArgumentParser(description='Remove water from within fiber cylinders')
    parser.add_argument('gro_input', help='Input GROMACS .gro file')
    parser.add_argument('top_input', help='Input GROMACS .top file')
    parser.add_argument('pdb_output', help='Output PDB file')
    parser.add_argument('-r', '--radius', type=float, default=2.0,
                       help='Cylinder radius (default: 2.0 nm)')
    
    args = parser.parse_args()
    
    # Read structure with MDAnalysis
    print(f"Reading {args.gro_input}...")
    u = mda.Universe(args.gro_input)
    
    print(f"Reading {args.top_input}...")
    top_lines = read_top_file(args.top_input)
    
    # Define cylinders
    atoms_per_cylinder = 660 * 28
    n_cylinders = 3
    
    # Calculate cylinder properties
    cylinders = []
    for i in range(n_cylinders):
        start_idx = i * atoms_per_cylinder
        end_idx = (i + 1) * atoms_per_cylinder
        
        # Get atoms for this cylinder
        cylinder_atoms = u.atoms[start_idx:end_idx]
        
        com_x, com_y, z_min, z_max = calculate_cylinder_properties(cylinder_atoms)
        cylinders.append({
            'com_x': com_x,
            'com_y': com_y,
            'z_min': z_min,
            'z_max': z_max
        })
        
        print(f"Cylinder {i+1}: COM=({com_x:.3f}, {com_y:.3f}), "
              f"z-range=[{z_min:.3f}, {z_max:.3f}]")
    
    # Select water molecules (assuming water has name W)
    water_atoms = u.select_atoms("name W")
    initial_water_count = len(water_atoms)
    
    # Find water molecules to remove
    water_indices_to_remove = []
    
    for water in water_atoms:
        water_pos = water.position
        
        # Check if water is inside any cylinder
        inside_cylinder = False
        for cyl in cylinders:
            if is_water_inside_cylinder(water_pos, cyl['com_x'], cyl['com_y'], 
                                      cyl['z_min'], cyl['z_max'], args.radius):
                inside_cylinder = True
                break
        
        if inside_cylinder:
            water_indices_to_remove.append(water.index)
    
    # Create selection of atoms to keep
    all_indices = set(range(len(u.atoms)))
    indices_to_remove = set(water_indices_to_remove)
    indices_to_keep = sorted(list(all_indices - indices_to_remove))
    
    # Create new universe with only atoms to keep
    atoms_to_keep = u.atoms[indices_to_keep]
    
    # Count remaining water molecules
    final_water_count = len(atoms_to_keep.select_atoms("name W"))
    removed_water_count = initial_water_count - final_water_count
    
    print(f"\nWater molecules: initial={initial_water_count}, "
          f"removed={removed_water_count}, final={final_water_count}")
    
    # Update topology file
    print(f"\nUpdating {args.top_input}...")
    updated_top_lines = update_water_count_in_top(top_lines, final_water_count)
    
    # Create backup of original topology
    backup_top = args.top_input + '.bak'
    write_top_file(backup_top, top_lines)
    print(f"Backup saved as {backup_top}")
    
    # Write updated topology
    write_top_file(args.top_input, updated_top_lines)
    
    # Write output PDB file
    print(f"\nWriting {args.pdb_output}...")
    
    # Write the atoms to keep to the PDB file
    atoms_to_keep.write(args.pdb_output)
    
    print(f"\nDone! Removed {removed_water_count} water molecules.")
    print(f"Output written to {args.pdb_output}")
    print(f"Topology updated in {args.top_input}")


if __name__ == "__main__":
    main()
