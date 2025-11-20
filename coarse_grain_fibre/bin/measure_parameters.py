#!/usr/bin/env python3
import argparse
import sys
import re
import os
import subprocess
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

NATOMS=28

def get_molecules(universe):
    """
    Get coarse-grained molecules from MDAnalysis universe when bond information is missing.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The MDAnalysis universe containing the trajectory
    atoms_per_molecule : int
        Number of atoms per molecule (default: 28)
    
    Returns:
    --------
    molecules : list of AtomGroup
        List of AtomGroup objects representing each molecule
    """    
    # Get total number of atoms
    n_atoms = len(universe.atoms)
    
    # Calculate number of molecules
    n_molecules = n_atoms // NATOMS
    
    if n_atoms % NATOMS != 0:
        print(f"Warning: {n_atoms} atoms is not divisible by {NATOMS}. "
              f"Last {n_atoms % NATOMS} atoms will be ignored.")
    
    molecules = []
    
    # Create AtomGroup for each molecule
    for i in range(n_molecules):
        start_idx = i * NATOMS
        end_idx = start_idx + NATOMS
        
        # Select atoms for this molecule
        mol_atoms = universe.atoms[start_idx:end_idx]
        molecules.append(mol_atoms)
    
    return molecules

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
        molecules = get_molecules(universe)
        xy_coms = []

        for mol in molecules:
            sel_atoms = mol.select_atoms(f"resname {' '.join(resnames)}")
            if len(sel_atoms) == 0:
                continue
            com = sel_atoms.center_of_geometry()
            xy_coms.append(com[:2])  # Project to XY

        if not xy_coms:
            print(f"Warning: No COMs found in frame {ts.frame} (t = {ts.time:.1f} ps)")
            continue

        xy_coms = np.array(xy_coms)
        centroid = xy_coms.mean(axis=0)
        distances = np.linalg.norm(xy_coms - centroid, axis=1)
        radius = distances.mean() + 2.75/2 # CG legs contain 1 bead that is considered "core" in AA. The mean distance is 2.75 (from CG distribution in water)

        radii.append(radius)
        times.append(ts.time)

        if i % max(1, universe.trajectory.n_frames // 10) == 0:
            all_coms_for_plot.append(xy_coms)

        # Simple check for system centering
        center_x, center_y = box[0] / 2, box[1] / 2
        margin_x, margin_y = box[0] * 0.05, box[1] * 0.05

        if not (center_x - margin_x <= centroid[0] <= center_x + margin_x and
                center_y - margin_y <= centroid[1] <= center_y + margin_y):
            continue
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
    molecules = get_molecules(universe)

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
            center = mol.center_of_geometry()
    
            # Calculate distances from center
            distances = np.linalg.norm(mol.positions - center, axis=1)
            
            # Calculate unweighted Rg
            rg = np.sqrt(np.mean(distances**2))
                    
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
        angle = [9,7,3,5]
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
    for i, ts in enumerate(universe.trajectory)[::10]:
        # Store time
        times.append(ts.time)
        
        # Get fragments (molecules)
        molecules = get_molecules(universe)
        
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
    import shutil
    rsrc = "../structure_files/vdwradii_CG.dat"
    rdst = "./vdwradii.dat"


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
    dt = 100e3
    
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
        "-o", str(temp_xvg),
        "-dt", str(dt)
    ])
    
    # Run gmx sasa with echo to provide the index group
    try:
        shutil.copy(rsrc, rdst)
        print(f"Copied {rsrc} → {rdst}")
        # Create input for group selection
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
        os.remove(rdst)
        print(f"Removed temporary file {rdst}")
        
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

def g_r(universe, resnames="MM2", r_cutoff=1.75, r_bin=0.2, last_quarter=True):
    """
    Compute g(r) for all molecule-molecule interactions on a cylindrical surface.

    Fixes:
    - Per-frame normalization: accumulate OBS and EXP per bin, then ratio.
    - Units are nm everywhere; output column is 'r (nm)'.
    - Safe geometric window enforced per frame: s_max = min(r_cutoff, 0.45*min(Lz, π*R)).
      Counts and expected beyond s_max are excluded that frame, preventing wrap artefacts.
    - Plot saved to ../logs/g_r/g_r.png
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import os

    os.makedirs("../logs/g_r", exist_ok=True)

    # Coerce inputs to float
    r_bin = float(r_bin)
    r_cutoff = float(r_cutoff)

    # Convert resnames to list if it's a string (allow comma-separated)
    if isinstance(resnames, str):
        resnames = [s.strip() for s in resnames.split(',') if s.strip()]

    # Determine frames to process
    n_frames = universe.trajectory.n_frames
    start_frame = int(0.75 * n_frames) if last_quarter else 0
    print(f"Processing frames {start_frame} to {n_frames - 1} out of {n_frames} total frames")

    # Resid IDs to analyze
    mol_resids = sorted(list(set(universe.residues.resids)))
    print(f"Molecule resids to analyze: {mol_resids}")

    # r bins (nm)
    r_bins = np.arange(0.0, r_cutoff + r_bin, r_bin)
    r_centers = r_bins[:-1] + r_bin / 2.0
    n_bins = len(r_centers)

    # Initialize pair bookkeeping
    interaction_pairs = []
    for i in mol_resids:
        for j in mol_resids:
            if i <= j:
                interaction_pairs.append((i, j))

    # Accumulators for counts (observed) and expected per bin
    counts = {f"{i}-{j}": np.zeros(n_bins) for (i, j) in interaction_pairs}
    expected = {f"{i}-{j}": np.zeros(n_bins) for (i, j) in interaction_pairs}

    total_processed_frames = 0
    radii_seen = []
    Lz_seen = []

    # You must provide get_molecules(universe) elsewhere (kept identical to your codebase)
    for ts in universe.trajectory[start_frame::10]:
        # Gather fragments and COMs per resid for the selected atoms
        fragments = get_molecules(universe)

        molecules_by_resid = {resid: [] for resid in mol_resids}
        coms_by_resid = {resid: [] for resid in mol_resids}
        temp_coms = []

        for frag in fragments:
            if len(frag.residues) == 0:
                continue
            resid = frag.residues[0].resid

            sel_atoms = frag.select_atoms(f"resname {' '.join(resnames)}")
            if len(sel_atoms) == 0:
                continue

            com = sel_atoms.center_of_geometry()
            temp_coms.append((resid, frag, com))

        if not temp_coms:
            print(f"Warning: No valid molecules found in frame {ts.frame}")
            continue

        # Center in xy-plane
        all_coms_array = np.array([com for _, _, com in temp_coms])
        xy_centroid = np.mean(all_coms_array[:, :2], axis=0)
        shift = np.array([-xy_centroid[0], -xy_centroid[1], 0.0])

        for resid, frag, com in temp_coms:
            centered = com + shift
            molecules_by_resid[resid].append(frag)
            coms_by_resid[resid].append(centered)

        # Compute radius (mean distance to z-axis) and box length along z
        all_coms = np.array([c for lst in coms_by_resid.values() for c in lst])
        if all_coms.size == 0:
            print(f"Warning: No molecules after centering in frame {ts.frame}")
            continue

        radius = float(np.mean(np.sqrt(all_coms[:, 0] ** 2 + all_coms[:, 1] ** 2)))
        box_z = float(universe.dimensions[2])

        radii_seen.append(radius)
        Lz_seen.append(box_z)

        # Geodesic distance on cylinder (θ,z) metric
        def com_distance(com_i, com_j, r, boxz):
            theta_i = np.arctan2(com_i[1], com_i[0])
            theta_j = np.arctan2(com_j[1], com_j[0])
            dtheta = np.abs(theta_i - theta_j)
            if dtheta > np.pi:
                dtheta = 2.0 * np.pi - dtheta
            arc = r * dtheta
            dz = np.abs(com_i[2] - com_j[2])
            if dz > boxz / 2.0:
                dz = boxz - dz
            return float(np.sqrt(arc * arc + dz * dz))

        # Per-frame densities
        frame_surface_area = 2.0 * np.pi * radius * box_z
        N_frame = {resid: len(coms_by_resid[resid]) for resid in mol_resids}
        rho_frame = {resid: (N_frame[resid] / frame_surface_area) if frame_surface_area > 0 else 0.0
                     for resid in mol_resids}

        # Locally flat validity window (per frame)
        s_flat = 0.45 * min(box_z, np.pi * radius)  # safety factor
        s_max = min(r_cutoff, s_flat)

        # Precompute valid bins mask for expected counts
        valid_bins = (r_centers <= s_max)
        shell_area = 2.0 * np.pi * r_centers * r_bin  # locally flat annulus area

        # Fill observed counts per pair
        for (i, j) in interaction_pairs:
            if not coms_by_resid[i] or not coms_by_resid[j]:
                continue

            coms_i = np.array(coms_by_resid[i])
            coms_j = np.array(coms_by_resid[j])

            key = f"{i}-{j}"

            for com_i in coms_i:
                for com_j in coms_j:
                    # Skip self for identical sets
                    if i == j and np.allclose(com_i, com_j):
                        continue
                    s = com_distance(com_i, com_j, radius, box_z)
                    if s > s_max:
                        continue
                    k = int(s / r_bin)
                    if 0 <= k < n_bins:
                        counts[key][k] += 1.0

            # Expected neighbors for each central i molecule in this frame (only within valid bins)
            if N_frame[i] > 0:
                exp_per_central = rho_frame[j] * shell_area  # array over bins
                # Zero out invalid bins for this frame
                exp_add = np.where(valid_bins, exp_per_central, 0.0)
                expected[key] += N_frame[i] * exp_add

        total_processed_frames += 1

    if total_processed_frames == 0:
        raise ValueError("No frames were processed. Check your data and filters.")

    # Final g(r): ratio of totals
    g_r_values = {}
    for key in counts.keys():
        denom = expected[key]
        mask = denom > 0.0
        g = np.zeros_like(denom)
        g[mask] = counts[key][mask] / denom[mask]
        g_r_values[key] = g

    # Create DataFrame (nm)
    df = pd.DataFrame({'r (Å)': r_centers})
    for pair, values in g_r_values.items():
        df[f'g(r) {pair}'] = values

    # Plot
    plt.figure(figsize=(10, 6))
    for pair in g_r_values.keys():
        plt.plot(r_centers, g_r_values[pair], label=f'{pair}')
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.title('Radial Distribution Functions on Cylindrical Surface (per-frame normalized)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("../logs/g_r/g_r.png", dpi=300)
    plt.close()

    # Report useful ranges
    if radii_seen and Lz_seen:
        r_mean = float(np.mean(radii_seen))
        Lz_min = float(np.min(Lz_seen))
        s_flat_global = 0.45 * min(Lz_min, np.pi * r_mean)
        print(f"Average cylinder radius ~ {r_mean:.3f} nm, min Lz ~ {Lz_min:.3f} nm.")
        print(f"Effective locally-flat window (global estimate) ≈ {s_flat_global:.3f} nm; "
              f"you set r_cutoff = {r_cutoff:.3f} nm.")

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
    'sasa' : sasa, 
    'g-r': g_r
    # Add more analysis methods here later
}
# Add this near the top of your file, before the analysis functions
def get_default_args(universe, structure_path):
    """Generate default arguments for analysis methods based on the universe data."""
    
    # Get all unique resids and exclude the highest one
    all_resids = sorted(list(set(universe.residues.resids)))
    default_resids_str = ','.join(str(r) for r in all_resids)

    # Get index file path by replacing extension of structure file
    structure_path_obj = Path(structure_path)
    index_file = str(structure_path_obj.with_suffix('.ndx'))

    return {
        'radius': {'resnames': "RLG,LLG"},
        'box-z' : {},
        'RG': {'resids': default_resids_str},
        'torsion': {
            'angles': '9,7,3,5', 
            'angle_names': 'theta'
        },  
        'sasa': {'index_file': index_file, 'index_group': 'MOTOR'},  
        'g-r': {'resnames': "MM2", 'r_cutoff': 20.0, 'r_bin': 0.2, 'last_quarter': True}
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
    universe = mda.Universe(args.structure, args.trajectory, topology_format='PDB')

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
        if len(func_kwargs) == 0:
            func_kwargs = default_args[args.analysis_method]
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