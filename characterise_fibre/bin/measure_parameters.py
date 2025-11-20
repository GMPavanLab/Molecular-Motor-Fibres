#!/usr/bin/env python3
import argparse
import sys
import re
import os
import subprocess
import shutil
from pathlib import Path
import time
import warnings
import pandas as pd
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import calc_dihedrals
import seaborn as sns
warnings.filterwarnings("ignore", message="Found no information for attr: 'formalcharges'")
warnings.filterwarnings("ignore", message="Found missing chainIDs.")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze molecular dynamics trajectories.')
    parser.add_argument('structure', help='Structure file (PDB)')
    parser.add_argument('trajectory', help='Trajectory file (XTC)')
    parser.add_argument('output_name', help='Base name for output files')
    parser.add_argument('analysis_method', help='Analysis method to apply or "ALL" to run all methods')
    parser.add_argument('--args', nargs='*', default=[], 
                        help='Optional arguments for the analysis method')
    
    return parser.parse_args()

def radius(universe, resnames=None):
    if resnames is None:
        raise ValueError("You must provide one or more RESNAMES via --args resnames=R1,R2,...")

    os.makedirs("../logs/radius", exist_ok=True)

    # Handle comma-separated list from string if needed
    if isinstance(resnames, str):
        resnames = resnames.split(',')

    # Get box dimensions from first frame for plotting
    universe.trajectory[0]
    box = universe.dimensions[:2]

    radii = []
    times = []
    all_coms_for_plot = []

    for i, ts in enumerate(universe.trajectory):
        molecules = universe.atoms.fragments
        xy_coms = []

        for mol in molecules:
            sel_atoms = mol.select_atoms(f"resname {' '.join(resnames)}")
            if len(sel_atoms) == 0:
                continue
            com = sel_atoms.center_of_mass()
            xy_coms.append(com[:2])  # Project to XY

        if not xy_coms:
            print(f"Warning: No COMs found in frame {ts.frame} (t = {ts.time:.1f} ps)")
            continue

        xy_coms = np.array(xy_coms)
        centroid = xy_coms.mean(axis=0)
        distances = np.linalg.norm(xy_coms - centroid, axis=1)
        radius = distances.mean()

        radii.append(radius)
        times.append(ts.time)

        if i % max(1, universe.trajectory.n_frames // 10) == 0:
            all_coms_for_plot.append(xy_coms)

        # Simple check for system centering
        center_x, center_y = box[0] / 2, box[1] / 2
        margin_x, margin_y = box[0] * 0.05, box[1] * 0.05

        if not (center_x - margin_x <= centroid[0] <= center_x + margin_x and
                center_y - margin_y <= centroid[1] <= center_y + margin_y):
            raise AssertionError(
                f"System appears not to be centered. Centroid at ({centroid[0]:.2f}, {centroid[1]:.2f}) "
                f"is outside the central 10% of the box [{center_x - margin_x:.2f}, {center_x + margin_x:.2f}]- "
                f"[{center_y - margin_y:.2f}, {center_y + margin_y:.2f}]"
            )

    # Save to CSV
    df = pd.DataFrame({'Time (ps)': times, 'Radius (Å)': radii})

    # Plot
    all_coms = np.vstack(all_coms_for_plot)
    plt.figure(figsize=(6, 6))
    sns.kdeplot(x=all_coms[:, 0], y=all_coms[:, 1], fill=True, cmap="viridis", bw_adjust=0.3, levels=100, thresh=0.05)

    # Draw box
    plt.plot([0, box[0], box[0], 0, 0], [0, 0, box[1], box[1], 0], color='red', lw=1)

    # Draw 1 nm reference scale in center
    mid = np.array(box) / 2
    plt.plot([mid[0] - 0.5, mid[0] + 0.5], [mid[1], mid[1]], color='white', lw=2)
    plt.text(mid[0], mid[1] + 0.2, '1 nm', ha='center', color='white')

    plt.xlabel("X (nm)")
    plt.ylabel("Y (nm)")
    plt.title("XY COM Density")
    plt.axis('equal')
    plt.tight_layout()
    base_name = Path(universe.filename).stem

    plot_path = f"../logs/radius/radius_check-{base_name}.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()

    return df

def RG(universe, resids=None):
    """
    Calculate radius of gyration for molecules of specified residue IDs.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The molecular system
    resids : str or None
        Comma-separated list of residue IDs to analyze
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing Rg values for each molecule at each timestep
    """
    if resids is None:
        raise ValueError("You must provide one or more RESIDS via --args resids=1,2,...")

    # Handle comma-separated list from string if needed
    if isinstance(resids, str):
        resids = [int(r) for r in resids.split(',')]
    
    # Create log directory if it doesn't exist
    os.makedirs("../logs/rg", exist_ok=True)
    
    # Dictionary to store Rg values for each fragment over time
    fragment_data = {}
    
    # List to store all Rg values for plotting
    all_rg_values = []
    resid_to_rg = {}  # For plotting: map resid to list of Rg values
    
    # Get relevant molecules (fragments) in the first frame to establish identities
    universe.trajectory[0]
    molecules = universe.atoms.fragments
    
    # Identify molecules that belong to specified resids
    target_molecules = []
    target_resids = []
    
    for mol in molecules:
        # Get the resid of the first atom in the molecule (all atoms in a fragment have the same resid)
        mol_resid = mol.residues[0].resid
        
        if mol_resid in resids:
            target_molecules.append(mol)
            target_resids.append(mol_resid)
            # Initialize entry in fragment_data dictionary
            fragment_data[mol.atoms[0].ix] = {'resid': mol_resid}
    
    if not target_molecules:
        raise ValueError(f"No molecules found with resids {resids}")
    
    # Collect time points to use as column names
    times = []
    
    # Track which frames to use for plotting (10% of frames)
    frames_for_plot = set(range(0, universe.trajectory.n_frames, 
                              max(1, universe.trajectory.n_frames // 10)))
    
    # Iterate through trajectory frames
    for i, ts in enumerate(universe.trajectory):
        times.append(ts.time)
        
        # Calculate Rg for each target molecule
        for mol, mol_resid in zip(target_molecules, target_resids):
            rg = mol.radius_of_gyration()
            
            # Store Rg value in fragment_data
            time_key = f"{ts.frame}"
            fragment_data[mol.atoms[0].ix][time_key] = rg       

            # Store data for plotting if this is a frame we want to include
            if i in frames_for_plot:
                all_rg_values.append((mol_resid, rg))
                
                if mol_resid not in resid_to_rg:
                    resid_to_rg[mol_resid] = []
                resid_to_rg[mol_resid].append(rg)
    
    # Create DataFrame from fragment_data
    df_data = []
    for frag_id, data in fragment_data.items():
        df_data.append(data)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(df_data)
        
    # Generate plot showing distribution of Rg for each resid
    plt.figure(figsize=(10, 6))
    
    # Convert plotting data to DataFrame for easier seaborn plotting
    plot_data = []
    for resid, rg_values in resid_to_rg.items():
        for rg in rg_values:
            plot_data.append({'ResID': f'ResID {resid}', 'Rg (Å)': rg})
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create density plot using seaborn
    sns.kdeplot(data=plot_df, x='Rg (Å)', hue='ResID', fill=True, common_norm=False, alpha=0.5)
    
    plt.title('Distribution of Radius of Gyration by Residue ID')
    plt.xlabel('Radius of Gyration (Å)')
    plt.ylabel('Density')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    base_name = Path(universe.filename).stem
    plt.savefig(f"../logs/rg/rg_check-{base_name}.png", dpi=300)
    plt.close()
    
    return results_df

def torsion(universe, angle=None, angle_name=None):
    """
    Calculate torsion angles for specified atoms in molecules.
    
    Parameters
    ----------
    universe : MDAnalysis Universe
        The loaded trajectory and topology
    angle : list or str, optional
        List of 4 atom indices to use for torsion calculation, relative to each molecule
        Default: [54, 64, 63, 109]
    angle_name : str, optional
        Name of the torsion angle for output files
        Default: None (will use "torsion")
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with torsion angles for each molecule at each timepoint
    """    
    # Set default torsion indices if not provided
    if angle is None:
        angle = [54,63,64,109]
    elif isinstance(angle, str):
        # Parse comma-separated string into a list of integers
        angle = [int(idx) for idx in angle.split(',')]
    
    # Set default angle name if not provided
    if angle_name is None:
        angle_name = "torsion"

    # Ensure we have exactly 4 atoms for the torsion calculation
    if len(angle) != 4:
        raise ValueError("Torsion calculation requires exactly 4 atom indices")
    
    os.makedirs("../logs/torsion", exist_ok=True)

    # Initialize data structures
    torsion_data = {}  # Dictionary to store torsion values by molecule ID
    times = []  # List to store frame timestamps
    
    # List to store torsion angles for plotting (10% of frames)
    plot_data = []
    
    # Process each frame
    for i, ts in enumerate(universe.trajectory):
        # Store time
        times.append(ts.time)
        
        # Get fragments (molecules)
        molecules = universe.atoms.fragments
        
        # Process each molecule that belongs to RESIDs 1-4
        for mol in molecules:
            # Check if molecule belongs to RESIDs 1-4
            if any(residue.resname == "NA" for residue in mol.residues):
                continue
                
            # Get the resid for this molecule
            mol_resid = mol.residues[0].resid  # Assuming all atoms in a fragment belong to same resid
            
            # Skip if number of atoms in molecule is less than max index in torsion
            if len(mol) <= max(angle):
                print(f"Warning: Molecule with RESID {mol_resid} has fewer atoms than required for torsion calculation")
                continue
            
            # Get the specified atoms for torsion calculation (using relative indices)
            torsion_atoms = [mol[idx] for idx in angle]
            
            # Calculate the torsion angle (in degrees)
            positions = np.array([atom.position for atom in torsion_atoms])
            angle_rad = calc_dihedrals(positions[0], positions[1], positions[2], positions[3], box=ts.dimensions)
            angle_deg = np.degrees(angle_rad)  # Convert to degrees
            # Check if it's an array and extract the value if needed
            if hasattr(angle_deg, '__len__'):
                angle_deg = angle_deg[0]
            
            # Store the torsion angle
            mol_key = f"{mol_resid}_{mol.ix[0]}"  # Use resid and first atom index as unique key
            if mol_key not in torsion_data:
                torsion_data[mol_key] = {"resid": mol_resid, "angles": []}
            
            torsion_data[mol_key]["angles"].append(angle_deg)
            
            # Store data for plotting (10% of frames)
            if i % max(1, universe.trajectory.n_frames // 10) == 0:
                plot_data.append({
                    "resid": mol_resid,
                    "angle": angle_deg,
                    "time": ts.time
                })
    
    # Create DataFrame with wide format (time as columns)
    df_data = []
    for mol_key, data in torsion_data.items():
        row = {"resid": data["resid"]}
        for i, t in enumerate(times):
            if i < len(data["angles"]):
                row[f"t{i+1}"] = data["angles"][i]
            else:
                row[f"t{i+1}"] = None
        df_data.append(row)
    
    df = pd.DataFrame(df_data)
    
    # Create the plot with distributions by RESID
    plot_df = pd.DataFrame(plot_data)
    
    if not plot_df.empty:
        plt.figure(figsize=(10, 6))

        # Shift angles from [-180, 180] to [0, 360] range
        for resid in sorted(plot_df["resid"].unique()):
            resid_data = plot_df[plot_df["resid"] == resid].copy()
            # Convert angles to [0, 360] range
            resid_data["angle_shifted"] = resid_data["angle"].apply(lambda x: x + 360 if x < 0 else x)
            
            sns.kdeplot(
                resid_data["angle_shifted"], 
                label=f"RESID {resid}",
                fill=True,
                alpha=0.3,
                common_norm=False
            )

        plt.xlabel("Torsion Angle (degrees)")
        plt.ylabel("Density")
        plt.title(f"Distribution of Torsion Angles {angle_name} for Atoms {angle}")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)

        # Add periodic boundary markers at 0 and 360
        plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
        plt.axvline(360, color='gray', linestyle='--', alpha=0.5)
        # Add a marker at 180 to show the center
        plt.axvline(180, color='red', linestyle='--', alpha=0.5)

        # Set x-axis to show the full range of torsion angles centered at 180
        plt.xlim(0, 360)
        plt.xticks([0, 60, 120, 180, 240, 300, 360], 
                ['0°', '60°', '120°', '180°', '240°', '300°', '360°'])

        # Save the plot
        base_name = Path(universe.filename).stem
        plot_path = f"../logs/torsion/torsion_{angle_name}-{base_name}.png"
        plt.savefig(plot_path, dpi=300)
        plt.close()
    
    return df

def voronota(universe, resids=None, resnames=None):
    """
    Calculate Voronota volumes and areas for molecules of specified residue IDs.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The molecular system
    resids : str or None
        Comma-separated list of residue IDs to analyze
    resnames : str or None
        Comma-separated list of residue names to filter the results
        
    Returns
    -------
    tuple of pd.DataFrame
        Two DataFrames containing volume and area values for each molecule at each timestep
    """   
    if resids is None:
        raise ValueError("You must provide one or more RESIDS via --args resids=1,2,...")

    # Handle comma-separated list from string if needed
    if isinstance(resids, str):
        resids = [int(r) for r in resids.split(',')]
    if isinstance(resnames, str):
        resnames = resnames.split(',')
    
    # Create log and temp directories if they don't exist
    os.makedirs("../logs/voronota", exist_ok=True)
    os.makedirs("../tmp", exist_ok=True)

    # Dictionaries to store volume and area values for each fragment over time
    volume_data = {}
    area_data = {}
    
    # Get relevant molecules (fragments) in the first frame to establish identities
    universe.trajectory[0]
    molecules = universe.atoms.fragments
    resid_filter = ",".join(str(r) for r in resids)
    
    # Identify molecules that belong to specified resids
    target_molecules = []
    target_resids = []
    molecule_indices = {}  # Map molecule to its index for consistent identification
    
    atom_to_mol = {}
    for mol_idx, mol in enumerate(molecules):
        mol_resid = mol.residues[0].resid
        if mol_resid in resids:
            target_molecules.append(mol)
            target_resids.append(mol_resid)
            molecule_indices[mol.atoms[0].ix+1] = mol_idx
            
            # Initialize entries in data dictionaries
            volume_data[mol_idx] = {'resid': mol_resid}
            area_data[mol_idx] = {'resid': mol_resid}
            
        # Add all atoms to the mapping regardless of resid (for contact calculations)
        for atom in mol:
            atom_to_mol[atom.id] = mol_idx
    if not target_molecules:
        raise ValueError(f"No molecules found with resids {resids}")
    
    # Collect time points to use as column names
    times = []
    
    # Track which frames to use for plotting (10% of frames)
    frames_for_plot = set(range(0, universe.trajectory.n_frames, 
                              max(1, universe.trajectory.n_frames // 10)))
    
    # Lists to store volume and area values for plotting
    plot_volumes = []
    plot_areas = []
    start_time = time.time()
    # Iterate through trajectory frames
    for i, ts in enumerate(universe.trajectory[::10]):
        if i > 0:
            elapsed = time.time() - start_time
            time_per_frame = elapsed / i
            frames_left = int(universe.trajectory.n_frames/10) - i
            eta = frames_left * time_per_frame
            eta_min = int(eta / 60)
            eta_sec = int(eta % 60)
            print(f"Processing frame {i+1}/{universe.trajectory.n_frames //10} - ETA: {eta_min}m {eta_sec}s")
        times.append(ts.time)
        
        # Create a temporary PDB file for the current frame
        temp_pdb = Path("../tmp/voronota_frame.pdb")
        universe.atoms.write(str(temp_pdb))
        
        try:
            # Run voronota and capture output
            voronota_vol_cmd = [
                "voronota-volumes", 
                "--input", str(temp_pdb), 
                "--input-filter-query", f"--match 'r<{resid_filter}>'"
            ]
            result = subprocess.run(voronota_vol_cmd, capture_output=True, text=True, check=True)
            voronota_output = result.stdout.splitlines()
            
            # Parse voronota output to get volumes per atom
            atom_volumes = {}
            for line in voronota_output:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        # Extract atom ID from the complex identifier
                        atom_id_str = parts[0]
                        # Find the atom number inside the a<NUMBER> part
                        atom_match = re.search(r'a<(\d+)>', atom_id_str)
                        if atom_match:
                            try:
                                atom_id = int(atom_match.group(1))
                                volume = float(parts[1])
                                atom_volumes[atom_id] = volume
                            except ValueError:
                                print(f"\nWarning: Parsing volume {atom_id_str} yielded ValueError")

                        else:
                            print(f"\nWarning: Could not parse {atom_id_str} in volume output")

            # Calculate contact areas between molecules using voronota-contacts
            # We want all atom-atom contacts, not just solvent-accessible ones
            voronota_contact_cmd = [
                "voronota-contacts", 
                "--input", str(temp_pdb), 
                "--input-filter-query", f"--match 'r<{resid_filter}>'"
            ]            
            
            area_result = subprocess.run(voronota_contact_cmd, capture_output=True, text=True, check=True)
            area_output = area_result.stdout.splitlines()
            
            # Parse to get areas per molecule
            mol_areas = {}
            for line in area_output:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 3:
                        # Extract atom IDs from complex identifiers
                        area = float(parts[2])  # Area is in the third column
                        atom1_match = re.search(r'a<(\d+)>', parts[0])
                        resname1_match = re.search(r'R<(\w+)>', parts[0])
                        atom2_match = re.search(r'a<(\d+)>', parts[1])
                        resname2_match = re.search(r'R<(\w+)>', parts[1])
                        mol1 = -1
                        mol2 = -1
                        atom1_valid_resname = False
                        atom2_valid_resname = False

                        try:  
                            if atom1_match and resname1_match:
                                atom1 = int(atom1_match.group(1))
                                resname1 = resname1_match.group(1)
                                mol1 = atom_to_mol.get(atom1)
                                atom1_valid_resname = resnames is None or resname1 in resnames
                            if atom2_match and resname2_match:
                                atom2 = int(atom2_match.group(1))
                                resname2 = resname2_match.group(1)
                                mol2 = atom_to_mol.get(atom2)
                                atom2_valid_resname = resnames is None or resname2 in resnames
                            if mol1 != mol2 or not (atom1_valid_resname and atom2_valid_resname):
                                if atom1_match and atom1_valid_resname:
                                    mol_areas[mol1] = mol_areas.get(mol1, 0) + area
                                if atom2_match and atom2_valid_resname:
                                    mol_areas[mol2] = mol_areas.get(mol2, 0) + area
                        except ValueError:
                            print(f"\nWarning: Parsing area {parts} yielded ValueError")
                    else:
                        print(f"\nWarning: Not enough elements in part {parts} in area output")

            # Now process volumes and areas for our target molecules
            for mol_idx in molecule_indices.values():
                mol = molecules[mol_idx]
                mol_resid = mol.residues[0].resid
                
                # Calculate total volume for this molecule
                if resnames is None:
                    # No filtering
                    total_volume = sum(atom_volumes.get(atom.id, 0) for atom in mol)
                else:
                    # Only include atoms with the specified resnames
                    total_volume = sum(atom_volumes.get(atom.id, 0) for atom in mol if atom.resname in resnames)
                
                # Get total contact area for this molecule
                total_area = mol_areas.get(mol_idx, 0)   
                
                # Store in data dictionaries
                time_key = f"t{i+1}"
                volume_data[mol_idx][time_key] = total_volume
                area_data[mol_idx][time_key] = total_area
                
                # Store for plotting if this is a selected frame
                if i in frames_for_plot:
                    plot_volumes.append({
                        'ResID': f'ResID {mol_resid}',
                        'Volume (Å³)': total_volume
                    })
                    plot_areas.append({
                        'ResID': f'ResID {mol_resid}',
                        'Area (Å²)': total_area
                    })
            
            
        except subprocess.CalledProcessError as e:
            print(f"Error running Voronota for frame {i}: {e}")
            print(f"Voronota stderr: {e.stderr}")
            continue
        finally:
            # Clean up temporary file
            if temp_pdb.exists():
                temp_pdb.unlink()
    print('Done')

    # Create DataFrames from collected data
    volume_df_data = []
    area_df_data = []
    
    for idx, data in volume_data.items():
        volume_df_data.append(data)
    
    for idx, data in area_data.items():
        area_df_data.append(data)
    
    volume_df = pd.DataFrame(volume_df_data)
    area_df = pd.DataFrame(area_df_data)
    
    # Create plots
    base_name = Path(universe.filename).stem

    if plot_volumes and plot_areas:
        # Volume distribution plot
        plt.figure(figsize=(10, 6))
        volume_plot_df = pd.DataFrame(plot_volumes)
        
        # Basic statistics for volumes to detect potential issues
        print("\nVolume statistics by ResID:")
        for resid in volume_plot_df['ResID'].unique():
            resid_data = volume_plot_df[volume_plot_df['ResID'] == resid]['Volume (Å³)']
            print(f"{resid}: mean={resid_data.mean():.2f}, std={resid_data.std():.2f}, min={resid_data.min():.2f}, max={resid_data.max():.2f}")
        
        sns.kdeplot(data=volume_plot_df, x='Volume (Å³)', hue='ResID', 
                   fill=True, common_norm=False, alpha=0.5)
        
        plt.title('Distribution of Molecular Volumes by Residue ID')
        plt.xlabel('Volume (Å³)')
        plt.ylabel('Density')
        plt.grid(alpha=0.3)
        plt.tight_layout()

        plt.savefig(f"../logs/voronota/volume_check-{base_name}.png", dpi=300)
        plt.close()
        
        # Area distribution plot
        plt.figure(figsize=(10, 6))
        area_plot_df = pd.DataFrame(plot_areas)
        
        # Basic statistics for areas to detect potential issues
        print("\nContact area statistics by ResID:")
        for resid in area_plot_df['ResID'].unique():
            resid_data = area_plot_df[area_plot_df['ResID'] == resid]['Area (Å²)']
            print(f"{resid}: mean={resid_data.mean():.2f}, std={resid_data.std():.2f}, min={resid_data.min():.2f}, max={resid_data.max():.2f}")
        
        sns.kdeplot(data=area_plot_df, x='Area (Å²)', hue='ResID', 
                   fill=True, common_norm=False, alpha=0.5)
        
        plt.title('Distribution of Molecular Contact Areas by Residue ID')
        plt.xlabel('Area (Å²)')
        plt.ylabel('Density')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(f"../logs/voronota/area_check-{base_name}.png", dpi=300)
        plt.close()
    
    return volume_df, area_df

def sasa(universe, index_file=None, index_group=None):
    """
    Calculate Solvent Accessible Surface Area (SASA) for each frame in the trajectory.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The molecular system
    index_file : str or None
        Path to GROMACS index file (.ndx)
    index_group : str or None
        Name of the index group to calculate SASA for
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing SASA values for each timestep
    """
    if index_file is None or index_group is None:
        raise ValueError("You must provide both index_file and index_group via --args index_file=path.ndx index_group=GroupName")
    
    # Create log directory if it doesn't exist
    os.makedirs("../logs/sasa", exist_ok=True)
    
    # Get the structure and trajectory paths from the universe
    structure_path = universe.filename
    structure_path_obj = Path(structure_path)
    index_file = str(structure_path_obj.with_suffix('.ndx'))
    trajectory_path = universe.trajectory.filename
    
    # Create temporary files for GROMACS output
    temp_xvg = "../tmp/sasa_temp.xvg"
    os.makedirs("../tmp", exist_ok=True)
    
    # Constant parameters for SASA calculation
    probe_size = 0.185  # nm
    ndots = 4800
    
    # Prepare the gmx sasa command
    sasa_cmd = [
        "gmx", "sasa",
        "-f", trajectory_path,
        "-s", structure_path
    ]
    
    # Add index file if provided
    if index_file:
        sasa_cmd.extend(["-n", index_file])
    
    # Add other parameters
    sasa_cmd.extend([
        "-ndots", str(ndots),
        "-probe", str(probe_size),
        "-o", temp_xvg
    ])
    vdw, dest = Path('../structure_files/vdwradii_AA.dat'), Path('vdwradii.dat')
    
    # Run gmx sasa with echo to provide the index group
    try:
        # Create input for group selection
        shutil.copy2(vdw, dest)
        shell_cmd = f"printf '{index_group}\\n' | " + " ".join(sasa_cmd)        
        # Run gmx sasa
        print(f"Running gmx sasa with index group: {index_group}")
        process = subprocess.run(
            shell_cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        dest.unlink(missing_ok=True)
        # Parse XVG file
        times = []
        sasa_values = []
        
        with open(temp_xvg, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and headers
                if line.startswith(('#', '@')):
                    continue
                    
                # Parse data lines
                parts = line.split()
                if len(parts) >= 2:
                    time = float(parts[0])  # Time in ps
                    sasa = float(parts[1])  # SASA in nm²
                    
                    times.append(time)
                    sasa_values.append(sasa)
        
        # Convert to DataFrame
        df = pd.DataFrame({'Time (ps)': times, 'SASA (nm²)': sasa_values})
        
        # Create visualization plot
        plt.figure(figsize=(10, 6))
        
        # Plot raw data
        plt.plot(df['Time (ps)'], df['SASA (nm²)'], alpha=0.5, label='Raw', color='lightblue')
        
        # Plot rolling average if we have enough frames
        if len(df) >= 5:
            rolling_avg = df['SASA (nm²)'].rolling(window=5, center=True).mean()
            rolling_std = df['SASA (nm²)'].rolling(window=5, center=True).std()
            
            plt.plot(df['Time (ps)'], rolling_avg, label='5-frame Rolling Avg', color='blue', linewidth=2)
            plt.fill_between(
                df['Time (ps)'], 
                rolling_avg - rolling_std,
                rolling_avg + rolling_std,
                alpha=0.3,
                color='blue',
                label='±1 Std Dev'
            )
        
        plt.xlabel('Time (ps)')
        plt.ylabel('SASA (nm²)')
        plt.title('Solvent Accessible Surface Area over Time')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        base_name = Path(universe.filename).stem
        # Save the plot
        plt.savefig(f"../logs/sasa/sasa_check-{base_name}.png", dpi=300)
        plt.close()
        
        return df
        
    except subprocess.CalledProcessError as e:
        print(f"Error running gmx sasa: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise
    finally:
        # Clean up temporary files
        if os.path.exists(temp_xvg):
            os.remove(temp_xvg)

def g_r(universe, resnames="COR", r_cutoff=20.0, r_bin=0.2, last_quarter=True):
    """
    Compute g(r) for all molecule-molecule interactions on a cylindrical surface.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The loaded MD universe containing the trajectory
    resnames : str or list, optional
        Residue name(s) to use for COM calculation
    r_cutoff : float, optional
        Maximum distance in nm for g(r) calculation
    r_bin : float, optional
        Bin size in nm for g(r) histogram
    last_quarter : bool, optional
        If True, only use the last 25% of frames
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with r values and g(r) for each interaction type
    """
    # Create output directory
    os.makedirs("../logs/g_r", exist_ok=True)
    r_bin=float(r_bin)
    r_cutoff=float(r_cutoff)
    # Convert resnames to list if it's a string
    if isinstance(resnames, str):
        resnames = resnames.split(',')
    
    # Determine number of frames to process
    n_frames = universe.trajectory.n_frames
    start_frame = int(0.75 * n_frames) if last_quarter else 0
    print(f"Processing frames {start_frame} to {n_frames-1} out of {n_frames} total frames")
    
    # Get the highest resid (NA) to exclude
    all_resids = sorted(list(set(universe.residues.resids)))
    na_resid = max(all_resids)
    print(f"Excluding NA ions with resid {na_resid}")
    
    # Get molecule resids (all except NA)
    mol_resids = [r for r in all_resids if r != na_resid]
    print(f"Molecule resids to analyze: {mol_resids}")
    
    # Set up the r bins for g(r)
    r_bins = np.arange(0, r_cutoff + r_bin, r_bin)
    r_centers = r_bins[:-1] + r_bin/2
    
    # Initialize dictionaries to store histograms and counts
    g_r_histograms = {}
    mol_counts = {resid: 0 for resid in mol_resids}
    
    def com_distance(com_i, com_j, r, box_z):
         # Project to cylinder surface
        theta_i = np.arctan2(com_i[1], com_i[0])
        z_i = com_i[2]
        theta_j = np.arctan2(com_j[1], com_j[0])
        z_j = com_j[2]
        
        # Calculate angular separation (minimum angle considering periodicity)
        d_theta = np.abs(theta_i - theta_j)
        if d_theta > np.pi:
            d_theta = 2 * np.pi - d_theta
        
        # Calculate arc length on cylinder surface
        arc_length = r * d_theta
        
        # Calculate z-distance considering periodic boundary
        d_z = np.abs(z_i - z_j)
        if d_z > box_z / 2:
            d_z = box_z - d_z
        
        # Total distance on cylinder surface (pythagorean)
        return np.sqrt(arc_length**2 + d_z**2)
    
    # Create all possible interaction pairs
    interaction_pairs = []
    for i in mol_resids:
        for j in mol_resids:
            # Only add the pair if i <= j (ensures we get 1-1, 1-2, 2-2, etc., but not 2-1)
            if i <= j:
                interaction_pairs.append((i, j))
                g_r_histograms[f"{i}-{j}"] = np.zeros(len(r_bins) - 1)
        
    total_processed_frames = 0
    
    # Process each frame
    for ts in universe.trajectory[start_frame:]:
        # Get all fragments (molecules)
        fragments = universe.atoms.fragments
        
        # Dictionary to store molecules by resid
        molecules_by_resid = {resid: [] for resid in mol_resids}
        coms_by_resid = {resid: [] for resid in mol_resids}
        temp_coms = []

        # Group molecules by resid and calculate their COMs
        for frag in fragments:
            # Get resid for this fragment
            if len(frag.residues) == 0:
                continue
                
            resid = frag.residues[0].resid
            
            # Skip NA molecules
            if resid == na_resid:
                continue
                
            # Select atoms with specified resname for COM calculation
            sel_atoms = frag.select_atoms(f"resname {' '.join(resnames)}")
            
            if len(sel_atoms) == 0:
                continue
                
            # Calculate COM
            com = sel_atoms.center_of_mass()
            
            # Store molecule and temporary COM
            temp_coms.append((resid, frag, com))

        # Calculate the centroid of all COMs in the xy-plane
        if temp_coms:
            all_coms_array = np.array([com for _, _, com in temp_coms])
            xy_centroid = np.mean(all_coms_array[:, :2], axis=0)
            
            # Shift vector (only in xy-plane, z remains unchanged)
            shift = np.array([-xy_centroid[0], -xy_centroid[1], 0])
            
            # Second pass: store the shifted COMs
            for resid, frag, com in temp_coms:
                # Apply the shift to center at (0,0) in xy-plane
                centered_com = com + shift
                
                # Store the molecule and its centered COM
                molecules_by_resid[resid].append(frag)
                coms_by_resid[resid].append(centered_com)
        else:
            print(f"Warning: No valid molecules found in frame {ts.frame}")
         
        # Count molecules of each type
        for resid in mol_resids:
            mol_counts[resid] += len(coms_by_resid[resid])
        
        # Calculate cylinder radius for this frame (average distance from z-axis)
        all_coms = []
        for resid in mol_resids:
            all_coms.extend(coms_by_resid[resid])
        
        all_coms = np.array(all_coms)
        # Add after calculating all_coms
        if len(all_coms) == 0:
            print(f"Warning: No molecules found in frame {ts.frame}")
            continue
            
        # Calculate distances from z-axis
        radius = np.mean(np.sqrt(all_coms[:, 0]**2 + all_coms[:, 1]**2))
        
        box_z = universe.dimensions[2]

        # Calculate g(r) for each interaction pair
        for i, j in interaction_pairs:
            if not coms_by_resid[i] or not coms_by_resid[j]:
                continue
                
            # Convert to numpy arrays for vectorized operations
            coms_i = np.array(coms_by_resid[i])
            coms_j = np.array(coms_by_resid[j])
            
            # For each COM in i, calculate distances to all COMs in j
            for com_i in coms_i:                
                for com_j in coms_j:
                    # Skip self-interactions
                    if i == j and np.allclose(com_i, com_j):
                        continue

                    distance = com_distance(com_i, com_j, radius, box_z)
                    # Skip distances beyond cutoff
                    if distance > r_cutoff:
                        continue
                    
                    # Add to histogram
                    bin_idx = int(distance / r_bin)
                    if bin_idx < len(g_r_histograms[f"{i}-{j}"]):
                        g_r_histograms[f"{i}-{j}"][bin_idx] += 1
        
        total_processed_frames += 1
    
    if total_processed_frames == 0:
        raise ValueError("No frames were processed. Check your data and filters.")
    
    # Normalize histograms to get g(r)
    g_r_values = {}
    
    for i, j in interaction_pairs:
        key = f"{i}-{j}"
        # Average number of molecules of type j
        n_j = mol_counts[j] / total_processed_frames        
        # Average density of molecules of type j
        avg_surface_area = 2 * np.pi * radius * box_z  # Using the last radius, could use average instead
        avg_rho_j = n_j / avg_surface_area
        
        # Calculate g(r)
        g_r_values[key] = np.zeros_like(g_r_histograms[key])
        
        for k in range(len(r_centers)):
            # Expected number in a shell = rho * (2πr * dr)
            shell_area = 2 * np.pi * r_centers[k] * r_bin 
            expected_count = avg_rho_j * shell_area
            
            # Only normalize if expected count is not zero
            if expected_count > 0:
                # Average count per molecule of type i
                avg_count = g_r_histograms[key][k] / mol_counts[i]
                g_r_values[key][k] = avg_count / expected_count
    
    # Create DataFrame
    df = pd.DataFrame({'r (Å)': r_centers})
    
    for pair, values in g_r_values.items():
        df[f'g(r) {pair}'] = values
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    for pair in g_r_values.keys():
        plt.plot(r_centers, g_r_values[pair], label=f'{pair}')
    
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.title('Radial Distribution Functions on Cylindrical Surface')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_path = f"../logs/g_r/g_r.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    return df

def box_z(universe):
    """
    Compute the z-dimension of the simulation box over time.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis Universe object containing the trajectory.
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns 'Time (ps)' and 'Z-dimension (Å)' 
    """
    os.makedirs("../logs/box_z", exist_ok=True)
    
    # Initialize lists to store data
    z_dimensions = []
    times = []
    
    # Iterate through trajectory frames
    for ts in universe.trajectory:
        # Get box dimensions
        box = universe.dimensions
        
        # Store z-dimension and time
        z_dimensions.append(box[2])  # Z-dimension is the third element
        times.append(ts.time)
    
    # Create DataFrame
    df = pd.DataFrame({'Time (ps)': times, 'Z-dimension (Å)': z_dimensions})
    
    # Plot the z-dimension over time
    plt.figure(figsize=(10, 6))
    plt.plot(times, z_dimensions)
    plt.xlabel('Time (ps)')
    plt.ylabel('Z-dimension (Å)')
    plt.title('Box Z-dimension over Time')
    plt.grid(True)
    plt.tight_layout()
    
    # Save plot
    base_name = Path(universe.filename).stem
    plot_path = f"../logs/box_z/box_z_dimension-{base_name}.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    return df

# Dictionary mapping analysis method names to functions
ANALYSIS_METHODS = {
    'radius': radius, 
    'box-z' : box_z, 
    'RG': RG,
    'torsion' : torsion,
    'voronota': voronota,
    'sasa' : sasa, 
    'g-r': g_r
    # Add more analysis methods here later
}
# Add this near the top of your file, before the analysis functions
def get_default_args(universe, structure_path):
    """Generate default arguments for analysis methods based on the universe data."""
    
    # Get all unique resids and exclude the highest one
    all_resids = sorted(list(set(universe.residues.resids)))
    default_resids = all_resids[:-1] if len(all_resids) > 1 else all_resids
    print(default_resids)
    default_resids_str = ','.join(str(r) for r in default_resids)

    # Get index file path by replacing extension of structure file
    structure_path_obj = Path(structure_path)
    index_file = str(structure_path_obj.with_suffix('.ndx'))

    return {
        'radius': {'resnames': "RLG,LLG"},
        'box-z' : {},
        'RG': {'resids': default_resids_str},
        'torsion': {
            'angles': '55,63,64,110;62,63,64,110;42,63,64,110', 
            'angle_names': 'theta;phi;alpha'
        },  
        'voronota': {'resids': default_resids_str, 'resnames': "COR"},
        'sasa': {'index_file': index_file, 'index_group': 'MOTOR'},  # Assuming "Protein" is a common group
        'g-r': {'resnames': "COR", 'r_cutoff': 20.0, 'r_bin': 0.2, 'last_quarter': True}
      }  

def main():
    args = parse_arguments()
    # Output directories
    data_dir = "../data"
    log_dir = "../logs/"
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    # Parse any additional arguments for the analysis method
    func_kwargs = {}
    for arg in args.args:
        if '=' in arg:
            key, value = arg.split('=', 1)
            func_kwargs[key] = value
    
    # Load the universe
    print(f"Loading trajectory: {args.structure} and {args.trajectory}")
    universe = mda.Universe(args.structure, args.trajectory, to_guess=['bonds', 'masses'], topology_format='PDB', vdwradii={'Na': 0.0})

    # Get default arguments based on the loaded universe
    default_args = get_default_args(universe, args.structure)

    # Handle "ALL" option
    if args.analysis_method.upper() == "ALL":
        print("Running all analysis methods with default arguments...")
        for method_name, method_func in ANALYSIS_METHODS.items():
            print(f"\nRunning analysis: {method_name}")
            try:
                # Use default arguments for this method
                default_kwargs = default_args.get(method_name, {})
                print(f"Using default arguments: {default_kwargs}")
                
                # Special handling for torsion if multiple angles are defined
                if method_name == 'torsion' and 'angles' in default_kwargs and 'angle_names' in default_kwargs:
                    angles = default_kwargs['angles'].split(';')
                    angle_names = default_kwargs['angle_names'].split(';')
                    
                    for angle, angle_name in zip(angles, angle_names):
                        print(f"Processing torsion angle: {angle_name}")
                        # Create a copy of default kwargs and update with current angle and name
                        kwargs_copy = default_kwargs.copy()
                        kwargs_copy['angle'] = angle
                        kwargs_copy['angle_name'] = angle_name
                        # Remove the multi-angle parameters to avoid confusion
                        if 'angles' in kwargs_copy: del kwargs_copy['angles']
                        if 'angle_names' in kwargs_copy: del kwargs_copy['angle_names']
                        
                        results = method_func(universe, **kwargs_copy)
                        # Save results
                        output_file = f"../data/{args.output_name}_torsion-{angle_name}.csv"
                        results.to_csv(output_file, index=False)
                        print(f"Results saved to {output_file}")
                else:
                    results = method_func(universe, **default_kwargs)
                    
                    # Handle results saving for voronota
                    if method_name == 'voronota':
                        vol_results, area_results = results
                        # Save both volume and area results
                        vol_output_file = f"../data/{args.output_name}_voronota-volume.csv"
                        vol_results.to_csv(vol_output_file, index=False)
                        print(f"Volume results saved to {vol_output_file}")
                        
                        area_output_file = f"../data/{args.output_name}_voronota-area.csv"
                        area_results.to_csv(area_output_file, index=False)
                        print(f"Area results saved to {area_output_file}")
                    else:
                        # Save results for other analysis methods
                        suffix = f'_{method_name}.csv'
                        output_file = f"../data/{args.output_name}{suffix}"
                        results.to_csv(output_file, index=False)
                        print(f"Results saved to {output_file}")
            except Exception as e:
                print(f"Error during {method_name} analysis: {e}")
                print("Continuing with next method...")
                continue
    else:
        # Single analysis method
        if args.analysis_method not in ANALYSIS_METHODS:
            print(f"Error: Unknown analysis method '{args.analysis_method}'")
            print(f"Available methods: {', '.join(ANALYSIS_METHODS.keys())}")
            sys.exit(1)
        
        # Run the selected analysis method
        analysis_func = ANALYSIS_METHODS[args.analysis_method]
        print(f"Running analysis: {args.analysis_method}")

        if not func_kwargs:
            func_kwargs = default_args[args.analysis_method]
            print(f"Using arguments for {args.analysis_method}: {func_kwargs}")

        # Special handling for torsion with multiple angles
        if args.analysis_method == 'torsion' and 'angles' in func_kwargs and 'angle_names' in func_kwargs:
            angles = func_kwargs['angles'].split(';')
            angle_names = func_kwargs['angle_names'].split(';')
            
            if len(angles) != len(angle_names):
                print("Error: The number of angles must match the number of angle names")
                sys.exit(1)
                
            for angle, angle_name in zip(angles, angle_names):
                print(f"Processing torsion angle: {angle_name}")
                # Create a copy of func_kwargs and update with current angle and name
                kwargs_copy = func_kwargs.copy()
                kwargs_copy['angle'] = angle
                kwargs_copy['angle_name'] = angle_name
                # Remove the multi-angle parameters to avoid confusion
                if 'angles' in kwargs_copy: del kwargs_copy['angles']
                if 'angle_names' in kwargs_copy: del kwargs_copy['angle_names']
                
                try:
                    results = analysis_func(universe, **kwargs_copy)
                    output_file = f"../data/{args.output_name}_torsion-{angle_name}.csv"
                    results.to_csv(output_file, index=False)
                    print(f"Results saved to {output_file}")
                except Exception as e:
                    print(f"Error processing torsion angle {angle_name}: {e}")
                    continue
        else:
            # Execute the regular analysis
            try:
                results = analysis_func(universe, **func_kwargs)

                if args.analysis_method == 'voronota':
                    vol_results, area_results = results
                    
                    # Save volume results
                    vol_output_file = f"../data/{args.output_name}_voronota-volume.csv"
                    vol_results.to_csv(vol_output_file, index=False)
                    print(f"Volume results saved to {vol_output_file}")
                    
                    # Save area results  
                    area_output_file = f"../data/{args.output_name}_voronota-area.csv"
                    area_results.to_csv(area_output_file, index=False)
                    print(f"Area results saved to {area_output_file}")
                else:
                    # Original code for other analysis methods
                    suffix = f'_{args.analysis_method}.csv'
                    output_file = f"../data/{args.output_name}{suffix}"
                    results.to_csv(output_file, index=False)
                    print(f"Results saved to {output_file}")
            
            except Exception as e:
                print(f"Error during analysis: {e}")
                sys.exit(1)

if __name__ == "__main__":
    main()