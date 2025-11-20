#!/usr/bin/env python3
"""
Calculate diffusion coefficients with drift correction using MDAnalysis
"""

import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, msd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argparse

def remove_pbc_effects(universe, unwrap=True):
    """Remove periodic boundary condition artifacts using proper unwrapping"""
    if unwrap:
        # Remove jumps from PBC
        ag = universe.atoms
        
        # Get first frame as reference
        universe.trajectory[0]
        pos_prev = ag.positions.copy()
        
        # Store unfolded positions
        unfolded_positions = []
        unfolded_positions.append(pos_prev.copy())
        
        # Process each frame
        for i, ts in enumerate(universe.trajectory):
            if i == 0:
                continue  # Skip first frame
                
            pos_current = ag.positions.copy()
            
            # Calculate displacement from previous frame
            displacement = pos_current - pos_prev
            
            # Check for wrapping (jumps)
            box = ts.dimensions[:3]
            if box is not None:
                # Correct for jumps in each dimension
                for dim in range(3):
                    # Find jumps (displacements larger than half box size)
                    jumps = np.where(displacement[:, dim] > box[dim]/2)[0]
                    displacement[jumps, dim] -= box[dim]
                    
                    jumps = np.where(displacement[:, dim] < -box[dim]/2)[0]
                    displacement[jumps, dim] += box[dim]
            
            # Update positions by adding corrected displacement
            pos_unfolded = unfolded_positions[-1] + displacement
            unfolded_positions.append(pos_unfolded.copy())
            pos_prev = pos_current.copy()
        
        # Replace trajectory positions with unfolded positions
        for i, ts in enumerate(universe.trajectory):
            ag.positions = unfolded_positions[i]
    
    return universe

def compute_system_drift(universe):
    """Compute average system drift between frames"""
    drifts = []
    first_frame = True
    prev_com = None
    
    for ts in universe.trajectory:
        # Calculate center of mass for entire system
        com = universe.atoms.center_of_mass()
        
        if not first_frame:
            drift = com - prev_com
            drifts.append(drift)
        else:
            first_frame = False
        
        prev_com = com.copy()
    
    return np.array(drifts)

def compute_fragment_deviations(universe, selection="all"):
    """Compute COG deviation for each fragment/molecule"""
    # Create residue groups (as proxy for molecules/fragments)
    fragment_size = 28

    # Get all atoms from your selection
    all_atoms = universe.select_atoms(selection)

    # Create fragments based on atom indices
    fragments = []
    for i in range(0, len(all_atoms), fragment_size):
        # Extract 28 atoms at a time as a fragment
        fragment = all_atoms[i:i+fragment_size]
        fragments.append(fragment)
    
    n_frames = len(universe.trajectory)
    n_fragments = len(fragments)
    print(n_fragments, n_frames)
    # Store initial positions
    initial_cogs = np.zeros((n_fragments, 3))
    deviations = np.zeros((n_frames, n_fragments, 3))
    
    # Get initial COGs
    universe.trajectory[0]
    for i, frag in enumerate(fragments):
        initial_cogs[i] = frag.center_of_geometry()
    
    # Calculate deviations for each frame
    for t, ts in enumerate(universe.trajectory):
        for i, frag in enumerate(fragments):
            current_cog = frag.center_of_geometry()
            deviations[t, i] = current_cog - initial_cogs[i]
            if t > 1 and abs(deviations[t, i][2] - deviations[t-1, i][2]) > 50:
                print("Warning, periodic jump")
    
    return deviations, fragments

def compute_msd_components(deviations, drift=None):
    """Calculate mean squared displacement in x, y, z directions"""
    n_frames = deviations.shape[0]
    msd_x = np.zeros(n_frames)
    msd_y = np.zeros(n_frames)
    msd_z = np.zeros(n_frames)
    msd_total = np.zeros(n_frames)
    
    for t in range(n_frames):
        # Calculate drift-corrected positions
        if drift is not None and t > 0:
            # Cumulative drift correction
            cumulative_drift = np.cumsum(drift[:t], axis=0)[-1]
            corrected_deviations = deviations[t] - cumulative_drift
        else:
            corrected_deviations = deviations[t]
        
        # Calculate squared deviations
        msd_x[t] = np.mean(corrected_deviations[:, 0]**2)
        msd_y[t] = np.mean(corrected_deviations[:, 1]**2)
        msd_z[t] = np.mean(corrected_deviations[:, 2]**2)
        msd_total[t] = np.mean(np.sum(corrected_deviations**2, axis=1))
    
    return msd_x, msd_y, msd_z, msd_total

def fit_diffusion_coefficient(time, msd, dimension=1, start_fit=0.2, end_fit=0.8):
    """Fit diffusion coefficient from MSD"""
    # Select fitting range
    start_idx = int(start_fit * len(time))
    end_idx = int(end_fit * len(time))
    
    # Perform linear fit
    slope, intercept, r_value, _, _ = stats.linregress(
        time[start_idx:end_idx], 
        msd[start_idx:end_idx]
    )
    
    # Calculate diffusion coefficient
    # D = slope / (2 * dimension)
    diffusion_coeff = slope / (2 * dimension)
    
    return diffusion_coeff, slope, r_value**2

def plot_msd_results(time, msd_x, msd_y, msd_z, fit_results=None):
    """Plot MSD components and fits"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    
    # Plot all components
    ax1.plot(time, msd_x, 'b-', label='MSD_x', alpha=0.7)
    ax1.plot(time, msd_y, 'g-', label='MSD_y', alpha=0.7)
    ax1.plot(time, msd_z, 'r-', label='MSD_z', alpha=0.7)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('MSD (Å²)')
    ax1.set_title('Mean Squared Displacement Components')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Zoom in on z-direction
    ax2.plot(time, msd_z, 'r-', label='MSD_z', alpha=0.7)
    
    if fit_results:
        D_z, slope, r2 = fit_results
        start_idx = int(0.2 * len(time))
        end_idx = int(0.8 * len(time))
        fit_line = slope * time + (msd_z[start_idx] - slope * time[start_idx])
        ax2.plot(time[start_idx:end_idx], fit_line[start_idx:end_idx], 
                'k--', linewidth=2, label=f'Fit: D_z = {D_z:.3e} Å²/ps\nR² = {r2:.3f}')
    
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('MSD_z (Å²)')
    ax2.set_title('Z-Direction Diffusion with Linear Fit')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('msd_diffusion_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a separate plot for the linear fit
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(time, msd_z, 'r-', label='MSD_z', alpha=0.7)
    
    if fit_results:
        D_z, slope, r2 = fit_results
        start_idx = int(0.2 * len(time))
        end_idx = int(0.8 * len(time))
        fit_line = slope * time + (msd_z[start_idx] - slope * time[start_idx])
        ax.plot(time[start_idx:end_idx], fit_line[start_idx:end_idx], 
                'k--', linewidth=2, label=f'Linear Fit')
        ax.set_title(f'Z-Direction Diffusion\nD_z = {D_z:.3e} Å²/ps\nR² = {r2:.3f}')
    
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('MSD_z (Å²)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('z_diffusion_fit.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Calculate diffusion coefficients with drift correction')
    parser.add_argument('topology', help='Topology file (e.g., .gro, .pdb)')
    parser.add_argument('trajectory', help='Trajectory file (e.g., .xtc, .trr)')
    parser.add_argument('--selection', default='all', help='Atom selection (default: all)')
    parser.add_argument('--unwrap', action='store_true', help='Unwrap trajectory')
    parser.add_argument('--output', default='diffusion_results', help='Output prefix')
    
    args = parser.parse_args()
    
    # Load universe
    u = mda.Universe(args.topology, args.trajectory)
    
    # Remove periodic boundary effects if requested
    if args.unwrap:
        u = remove_pbc_effects(u, unwrap=True)
    
    # Compute system drift
    print("Computing system drift...")
    drift = compute_system_drift(u)
    
    # Compute fragment deviations
    print("Computing fragment deviations...")
    deviations, fragments = compute_fragment_deviations(u, args.selection)
    
    # Calculate MSD components with drift correction
    print("Calculating MSD components...")
    msd_x, msd_y, msd_z, msd_total = compute_msd_components(deviations, drift)
    
    # Get time array
    time = np.array([ts.time for ts in u.trajectory])
    
    # Fit diffusion coefficient for z-direction
    print("Fitting diffusion coefficient...")
    D_z, slope, r2 = fit_diffusion_coefficient(time, msd_z, dimension=1)
    
    # Plot results
    plot_msd_results(time, msd_x, msd_y, msd_z, (D_z, slope, r2))
    
    # Save numerical results
    np.savetxt(f'{args.output}_msd_data.txt', 
               np.column_stack((time, msd_x, msd_y, msd_z, msd_total)),
               header='Time(ps) MSD_x MSD_y MSD_z MSD_total',
               fmt='%.6f')
    
    # Save drift data
    np.savetxt(f'{args.output}_drift.txt', drift,
               header='Drift_x Drift_y Drift_z',
               fmt='%.6f')
    
    # Print results
    print(f"\nResults:")
    print(f"Z-direction diffusion coefficient: {D_z:.3e} Å²/ps")
    print(f"Z-direction diffusion coefficient: {D_z*1e4:.3f} nm²/mu s")
    print(f"Linear fit R²: {r2:.3f}")
    print(f"Number of fragments analyzed: {len(fragments)}")
    
    # Create summary plot of drift
    plt.figure(figsize=(10, 6))
    plt.plot(time[1:], drift[:, 0], 'b-', label='Drift_x', alpha=0.7)
    plt.plot(time[1:], drift[:, 1], 'g-', label='Drift_y', alpha=0.7)
    plt.plot(time[1:], drift[:, 2], 'r-', label='Drift_z', alpha=0.7)
    plt.xlabel('Time (ps)')
    plt.ylabel('System Drift (Å)')
    plt.title('System Center of Mass Drift')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{args.output}_drift_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create additional analysis plots
    plt.figure(figsize=(12, 8))
    
    # Plot cumulative drift effect
    plt.subplot(2, 2, 1)
    cumulative_drift = np.cumsum(drift, axis=0)
    plt.plot(time[1:], cumulative_drift[:, 2], 'r-', label='Cumulative Z-drift')
    plt.xlabel('Time (ps)')
    plt.ylabel('Cumulative Drift (Å)')
    plt.title('Cumulative Z-direction Drift')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot MSD comparison
    plt.subplot(2, 2, 2)
    plt.plot(time, msd_z, 'r-', label='With drift correction', alpha=0.7)
    # Calculate uncorrected MSD for comparison
    msd_uncorrected = np.mean(deviations[:, :, 2]**2, axis=1)
    plt.plot(time, msd_uncorrected, 'b--', label='Without drift correction', alpha=0.7)
    plt.xlabel('Time (ps)')
    plt.ylabel('MSD_z (Å²)')
    plt.title('Drift Correction Effect on MSD_z')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot fragment distribution
    plt.subplot(2, 2, 3)
    final_positions = deviations[-1, :, 2]
    plt.hist(final_positions, bins=30, alpha=0.7, color='green')
    plt.xlabel('Z-displacement (Å)')
    plt.ylabel('Count')
    plt.title('Fragment Z-displacement Distribution')
    plt.grid(True, alpha=0.3)
    
    # Plot anisotropy
    plt.subplot(2, 2, 4)
    anisotropy = msd_z / ((msd_x + msd_y) / 2)
    plt.plot(time, anisotropy, 'purple', alpha=0.7)
    plt.axhline(y=1, color='k', linestyle='--', alpha=0.5)
    plt.xlabel('Time (ps)')
    plt.ylabel('Anisotropy (D_z/D_xy)')
    plt.title('Diffusion Anisotropy')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{args.output}_detailed_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()