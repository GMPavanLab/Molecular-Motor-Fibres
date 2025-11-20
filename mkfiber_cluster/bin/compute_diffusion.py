#!/usr/bin/env python3
"""
Calculate diffusion coefficients with best practices implementation
Based on LiveCoMS Best Practices Guide (Maginn et al., 2019)

Requirements:
- Input trajectory must be unwrapped using: gmx trjconv -mol -nojump
- For accurate results, provide system viscosity for finite size correction
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argparse
import warnings

# Physical constants
KB = 1.380649e-23  # Boltzmann constant (J/K)
XI = 2.837298      # Dimensionless constant for Yeh-Hummer correction

def validate_unwrapped_trajectory(universe, fragment_size=28):
    """
    Validate that input trajectory is properly unwrapped
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The trajectory universe
    fragment_size : int
        Expected size of molecular fragments
        
    Returns:
    --------
    bool : True if trajectory appears unwrapped
    """
    print("Validating trajectory unwrapping...")
    
    # Check for reasonable molecular displacements
    universe.trajectory[0]
    pos_first = universe.atoms.positions.copy()
    universe.trajectory[-1]
    pos_last = universe.atoms.positions.copy()
    
    max_displacement = np.max(np.linalg.norm(pos_last - pos_first, axis=1))
    box_size = np.min(universe.dimensions[:3])
    
    if max_displacement < box_size / 4:
        warnings.warn(
            f"Maximum displacement ({max_displacement:.2f} Å) seems small relative to "
            f"box size ({box_size:.2f} Å). Trajectory may still be wrapped.\n"
            "Please ensure trajectory is unwrapped using: gmx trjconv -mol -nojump"
        )
        return False
    
    print(f"✓ Maximum displacement: {max_displacement:.2f} Å")
    print(f"✓ Box size: {box_size:.2f} Å")
    return True

def validate_simulation_quality(universe, temperature_tolerance=5.0):
    """
    Perform quality checks on the simulation
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The trajectory universe
    temperature_tolerance : float
        Acceptable temperature deviation (K)
    """
    print("\nPerforming simulation quality checks...")
    
    n_frames = len(universe.trajectory)
    print(f"✓ Trajectory length: {n_frames} frames")
    
    # Check simulation box stability
    dimensions = []
    for ts in universe.trajectory:
        dimensions.append(ts.dimensions[:3])
    dimensions = np.array(dimensions)
    
    box_drift = np.std(dimensions, axis=0) / np.mean(dimensions, axis=0)
    if np.any(box_drift > 0.01):
        warnings.warn("Significant box size fluctuations detected (>1%)")
    else:
        print("✓ Box dimensions stable")
    
    return True

def compute_system_drift(universe):
    """
    Compute average system drift between frames
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The trajectory universe
        
    Returns:
    --------
    np.ndarray : Array of drift vectors [n_frames-1, 3]
    """
    print("Computing system center-of-mass drift...")
    
    drifts = []
    prev_com = None
    
    for ts in universe.trajectory:
        com = universe.atoms.center_of_mass()
        
        if prev_com is not None:
            drift = com - prev_com
            drifts.append(drift)
        
        prev_com = com.copy()
    
    drifts = np.array(drifts)
    
    # Report drift statistics
    drift_magnitude = np.linalg.norm(drifts, axis=1)
    print(f"✓ Average drift per frame: {np.mean(drift_magnitude):.4f} Å")
    print(f"✓ Maximum drift per frame: {np.max(drift_magnitude):.4f} Å")
    
    return drifts

def compute_fragment_deviations(universe, selection="all", fragment_size=28):
    """
    Compute center-of-geometry deviations for molecular fragments
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The trajectory universe
    selection : str
        Atom selection string
    fragment_size : int
        Number of atoms per fragment
        
    Returns:
    --------
    tuple : (deviations array [n_frames, n_fragments, 3], fragment list)
    """
    print(f"Computing fragment deviations (fragment size: {fragment_size})...")
    
    all_atoms = universe.select_atoms(selection)
    
    # Create fragments
    fragments = []
    for i in range(0, len(all_atoms), fragment_size):
        fragment = all_atoms[i:i+fragment_size]
        if len(fragment) == fragment_size:  # Only include complete fragments
            fragments.append(fragment)
    
    n_frames = len(universe.trajectory)
    n_fragments = len(fragments)
    
    print(f"✓ Analyzing {n_fragments} complete fragments")
    
    # Initialize arrays
    initial_cogs = np.zeros((n_fragments, 3))
    deviations = np.zeros((n_frames, n_fragments, 3))
    
    # Get initial center-of-geometry positions
    universe.trajectory[0]
    for i, frag in enumerate(fragments):
        initial_cogs[i] = frag.center_of_geometry()
    
    # Calculate deviations for each frame
    for t, ts in enumerate(universe.trajectory):
        for i, frag in enumerate(fragments):
            current_cog = frag.center_of_geometry()
            deviations[t, i] = current_cog - initial_cogs[i]
    
    return deviations, fragments

def validate_diffusive_regime(msd, time, box_length, estimated_radius=5.0):
    """
    Validate that simulation captures diffusive regime
    
    Parameters:
    -----------
    msd : np.ndarray
        Mean squared displacement data
    time : np.ndarray
        Time points
    box_length : float
        Simulation box length
    estimated_radius : float
        Estimated molecular radius (Å)
        
    Returns:
    --------
    bool : True if diffusive regime is captured
    """
    print("\nValidating diffusive regime...")
    
    final_rms = np.sqrt(msd[-1])
    
    # Check 1: RMS displacement > molecular radius
    if final_rms <= estimated_radius:
        warnings.warn(
            f"RMS displacement ({final_rms:.2f} Å) ≤ estimated molecular radius ({estimated_radius:.2f} Å). "
            "Simulation may be too short."
        )
        return False
    
    # Check 2: RMS displacement should be significant fraction of box
    if final_rms < box_length / 10:
        warnings.warn(
            f"RMS displacement ({final_rms:.2f} Å) < 10% of box length ({box_length:.2f} Å). "
            "Consider longer simulation for better statistics."
        )
    
    # Check 3: Verify linear regime in log-log plot
    log_time = np.log10(time[1:])  # Skip t=0
    log_msd = np.log10(msd[1:])
    
    # Fit slope in latter half
    mid_point = len(log_time) // 2
    slope, _, r_value, _, _ = stats.linregress(log_time[mid_point:], log_msd[mid_point:])
    
    print(f"✓ Final RMS displacement: {final_rms:.2f} Å")
    print(f"✓ Box length: {box_length:.2f} Å")
    print(f"✓ Log-log slope in diffusive regime: {slope:.3f} (should ≈ 1.0)")
    print(f"✓ R² for diffusive regime: {r_value**2:.3f}")
    
    if abs(slope - 1.0) > 0.2:
        warnings.warn(f"Log-log slope ({slope:.3f}) deviates significantly from 1.0")
        return False
    
    return True

def compute_msd_components(deviations, drift=None):
    """
    Calculate mean squared displacement with drift correction
    
    Parameters:
    -----------
    deviations : np.ndarray
        Fragment deviation array [n_frames, n_fragments, 3]
    drift : np.ndarray, optional
        System drift correction [n_frames-1, 3]
        
    Returns:
    --------
    tuple : (msd_x, msd_y, msd_z, msd_total)
    """
    n_frames = deviations.shape[0]
    msd_x = np.zeros(n_frames)
    msd_y = np.zeros(n_frames)
    msd_z = np.zeros(n_frames)
    msd_total = np.zeros(n_frames)
    
    for t in range(n_frames):
        # Apply drift correction
        corrected_deviations = deviations[t].copy()
        
        if drift is not None and t > 0:
            # Cumulative drift correction
            cumulative_drift = np.sum(drift[:t], axis=0)
            corrected_deviations -= cumulative_drift
        
        # Calculate MSD components
        msd_x[t] = np.mean(corrected_deviations[:, 0]**2)
        msd_y[t] = np.mean(corrected_deviations[:, 1]**2)
        msd_z[t] = np.mean(corrected_deviations[:, 2]**2)
        msd_total[t] = np.mean(np.sum(corrected_deviations**2, axis=1))
    
    return msd_x, msd_y, msd_z, msd_total

def plot_drift_correction_validation(time, deviations, drift, output_prefix="validation"):
    """
    Plot the effect of drift correction on molecular trajectories
    
    Parameters:
    -----------
    time : np.ndarray
        Time points
    deviations : np.ndarray
        Fragment deviation array [n_frames, n_fragments, 3]
    drift : np.ndarray
        System drift correction [n_frames-1, 3]
    output_prefix : str
        Output file prefix
    """
    print("Creating drift correction validation plots...")
    
    # Calculate MSD with and without drift correction
    n_frames = deviations.shape[0]
    
    # Without drift correction
    msd_y_uncorrected = np.zeros(n_frames)
    for t in range(n_frames):
        msd_y_uncorrected[t] = np.mean(deviations[t, :, 1]**2)  # Y-axis (index 1)
    
    # With drift correction
    msd_y_corrected = np.zeros(n_frames)
    for t in range(n_frames):
        corrected_deviations = deviations[t].copy()
        if t > 0:
            cumulative_drift = np.sum(drift[:t], axis=0)
            corrected_deviations -= cumulative_drift
        msd_y_corrected[t] = np.mean(corrected_deviations[:, 1]**2)
    
    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. MSD comparison
    ax1 = axes[0, 0]
    ax1.plot(time, msd_y_uncorrected, 'r-', label='Without drift correction', alpha=0.7)
    ax1.plot(time, msd_y_corrected, 'b-', label='With drift correction', alpha=0.7)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('MSD_y (Å²)')
    ax1.set_title('Effect of Drift Correction on MSD_y')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Cumulative drift vs time
    ax2 = axes[0, 1]
    cumulative_drift_y = np.zeros(n_frames)
    for t in range(1, n_frames):
        cumulative_drift_y[t] = np.sum(drift[:t, 1])  # Y-component
    
    ax2.plot(time, cumulative_drift_y, 'g-', linewidth=2)
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Cumulative Drift_y (Å)')
    ax2.set_title('System Drift Accumulation (Y-axis)')
    ax2.grid(True, alpha=0.3)
    
    # 3. Drift magnitude over time
    ax3 = axes[1, 0]
    drift_magnitude = np.linalg.norm(drift, axis=1)
    ax3.plot(time[1:], drift_magnitude, 'purple', alpha=0.7)
    ax3.axhline(y=np.mean(drift_magnitude), color='red', linestyle='--', 
               label=f'Mean: {np.mean(drift_magnitude):.3f} Å')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Drift Magnitude per Frame (Å)')
    ax3.set_title('System Drift Magnitude')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Ratio of drift to molecular displacement
    ax4 = axes[1, 1]
    # Calculate molecular displacement per frame (RMS of all molecules)
    molecular_displacement = np.zeros(n_frames-1)
    for t in range(1, n_frames):
        displacement = deviations[t] - deviations[t-1]
        molecular_displacement[t-1] = np.sqrt(np.mean(np.sum(displacement**2, axis=1)))
    
    drift_to_molecule_ratio = drift_magnitude / molecular_displacement
    ax4.plot(time[1:], drift_to_molecule_ratio, 'orange', alpha=0.7)
    ax4.axhline(y=np.mean(drift_to_molecule_ratio), color='red', linestyle='--',
               label=f'Mean ratio: {np.mean(drift_to_molecule_ratio):.1f}')
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Drift/Molecular Displacement Ratio')
    ax4.set_title('Drift vs Molecular Motion Ratio')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_drift_validation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print validation statistics
    print(f"\nDrift Correction Statistics:")
    print(f"  Average drift magnitude: {np.mean(drift_magnitude):.3f} Å/frame")
    print(f"  Average molecular displacement: {np.mean(molecular_displacement):.3f} Å/frame")
    print(f"  Drift/molecular ratio: {np.mean(drift_to_molecule_ratio):.1f}")
    print(f"  Final cumulative drift (Y): {cumulative_drift_y[-1]:.2f} Å")
    print(f"  MSD_y reduction due to drift correction: {(msd_y_uncorrected[-1] - msd_y_corrected[-1])/msd_y_uncorrected[-1]*100:.1f}%")

def plot_multiple_molecule_trajectories(deviations, drift, time, output_prefix="validation", molecule_indices=[0]):
    """
    Plot random walk of multiple molecules along Y-axis and in XY plane
    
    Parameters:
    -----------
    deviations : np.ndarray
        Fragment deviation array [n_frames, n_fragments, 3]
    drift : np.ndarray
        System drift correction [n_frames-1, 3]
    time : np.ndarray
        Time points
    output_prefix : str
        Output file prefix
    molecule_indices : list
        Indices of molecules to track
    """
    print(f"Creating multiple molecule trajectory plots for molecules: {molecule_indices}")
    
    n_frames = deviations.shape[0]
    n_molecules = len(molecule_indices)
    
    # Color palette for different molecules
    colors = plt.cm.tab10(np.linspace(0, 1, n_molecules))
    
    # Create plots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    
    # Storage for statistics
    all_steps_uncorrected = []
    all_steps_corrected = []
    molecule_stats = []
    
    for i, mol_idx in enumerate(molecule_indices):
        if mol_idx >= deviations.shape[1]:
            print(f"Warning: Molecule index {mol_idx} exceeds available molecules ({deviations.shape[1]})")
            continue
            
        color = colors[i]
        
        # Extract molecule trajectory
        mol_trajectory = deviations[:, mol_idx, :].copy()  # [n_frames, 3]
        
        # Apply drift correction
        mol_trajectory_corrected = mol_trajectory.copy()
        for t in range(1, n_frames):
            cumulative_drift = np.sum(drift[:t], axis=0)
            mol_trajectory_corrected[t] -= cumulative_drift
        
        # 1. Y-axis trajectory (uncorrected)
        ax1 = axes[0, 0]
        ax1.plot(time, mol_trajectory[:, 1], color=color, alpha=0.7, linewidth=1, 
                label=f'Mol {mol_idx}' if i < 5 else "")  # Limit legend entries
        
        # 2. Y-axis trajectory (drift corrected)
        ax2 = axes[0, 1]
        ax2.plot(time, mol_trajectory_corrected[:, 1], color=color, alpha=0.7, linewidth=1,
                label=f'Mol {mol_idx}' if i < 5 else "")
        
        # 3. Y-axis comparison for first molecule only (to avoid clutter)
        if i == 0:
            ax3 = axes[0, 2]
            ax3.plot(time, mol_trajectory[:, 1], 'r-', alpha=0.7, label='Uncorrected', linewidth=2)
            ax3.plot(time, mol_trajectory_corrected[:, 1], 'b-', alpha=0.7, label='Drift corrected', linewidth=2)
        
        # 4. XY trajectory (uncorrected)
        ax4 = axes[1, 0]
        ax4.plot(mol_trajectory[:, 0], mol_trajectory[:, 1], color=color, alpha=0.6, linewidth=1)
        ax4.scatter(mol_trajectory[0, 0], mol_trajectory[0, 1], color=color, s=30, 
                   marker='o', alpha=0.8, edgecolors='black', linewidths=0.5)
        ax4.scatter(mol_trajectory[-1, 0], mol_trajectory[-1, 1], color=color, s=30, 
                   marker='s', alpha=0.8, edgecolors='black', linewidths=0.5)
        
        # 5. XY trajectory (drift corrected)
        ax5 = axes[1, 1]
        ax5.plot(mol_trajectory_corrected[:, 0], mol_trajectory_corrected[:, 1], color=color, alpha=0.6, linewidth=1)
        ax5.scatter(mol_trajectory_corrected[0, 0], mol_trajectory_corrected[0, 1], color=color, s=30, 
                   marker='o', alpha=0.8, edgecolors='black', linewidths=0.5)
        ax5.scatter(mol_trajectory_corrected[-1, 0], mol_trajectory_corrected[-1, 1], color=color, s=30, 
                   marker='s', alpha=0.8, edgecolors='black', linewidths=0.5)
        
        # Calculate step sizes for histogram
        steps_uncorrected = np.diff(mol_trajectory[:, 1])  # Y-direction steps
        steps_corrected = np.diff(mol_trajectory_corrected[:, 1])
        all_steps_uncorrected.extend(steps_uncorrected)
        all_steps_corrected.extend(steps_corrected)
        
        # Calculate statistics for this molecule
        y_displacement_uncorrected = mol_trajectory[-1, 1] - mol_trajectory[0, 1]
        y_displacement_corrected = mol_trajectory_corrected[-1, 1] - mol_trajectory_corrected[0, 1]
        
        rms_displacement_uncorrected = np.sqrt(np.mean(mol_trajectory[-1]**2))
        rms_displacement_corrected = np.sqrt(np.mean(mol_trajectory_corrected[-1]**2))
        
        # Estimate diffusion coefficient for this molecule
        y_msd_single = mol_trajectory_corrected[:, 1]**2
        D_single = np.nan
        r_val_sq = np.nan
        if len(time) > 10:
            start_fit = int(0.2 * len(time))
            try:
                slope, _, r_val, _, _ = stats.linregress(time[start_fit:], y_msd_single[start_fit:])
                D_single = slope / 2  # 1D diffusion
                r_val_sq = r_val**2
            except:
                pass
        
        molecule_stats.append({
            'idx': mol_idx,
            'y_disp_uncorr': y_displacement_uncorrected,
            'y_disp_corr': y_displacement_corrected,
            'rms_uncorr': rms_displacement_uncorrected,
            'rms_corr': rms_displacement_corrected,
            'D_single': D_single,
            'r_squared': r_val_sq,
            'step_rms_uncorr': np.sqrt(np.mean(steps_uncorrected**2)),
            'step_rms_corr': np.sqrt(np.mean(steps_corrected**2))
        })
    
    # Set up axis labels and titles
    axes[0, 0].set_xlabel('Time (ps)')
    axes[0, 0].set_ylabel('Y displacement (Å)')
    axes[0, 0].set_title(f'Multiple Molecule Y-trajectories (Uncorrected)\n{n_molecules} molecules')
    axes[0, 0].grid(True, alpha=0.3)
    if n_molecules <= 5:
        axes[0, 0].legend()
    
    axes[0, 1].set_xlabel('Time (ps)')
    axes[0, 1].set_ylabel('Y displacement (Å)')
    axes[0, 1].set_title(f'Multiple Molecule Y-trajectories (Drift Corrected)\n{n_molecules} molecules')
    axes[0, 1].grid(True, alpha=0.3)
    if n_molecules <= 5:
        axes[0, 1].legend()
    
    axes[0, 2].set_xlabel('Time (ps)')
    axes[0, 2].set_ylabel('Y displacement (Å)')
    axes[0, 2].set_title(f'Y-trajectory Comparison\n(Molecule {molecule_indices[0]} only)')
    axes[0, 2].legend()
    axes[0, 2].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('X displacement (Å)')
    axes[1, 0].set_ylabel('Y displacement (Å)')
    axes[1, 0].set_title(f'XY Random Walks (Uncorrected)\nCircles=start, Squares=end')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axis('equal')
    
    axes[1, 1].set_xlabel('X displacement (Å)')
    axes[1, 1].set_ylabel('Y displacement (Å)')
    axes[1, 1].set_title(f'XY Random Walks (Drift Corrected)\nCircles=start, Squares=end')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axis('equal')
    
    # 6. Combined step size distribution
    ax6 = axes[1, 2]
    ax6.hist(all_steps_uncorrected, bins=50, alpha=0.5, color='red', 
             label=f'Uncorrected ({len(all_steps_uncorrected)} steps)', density=True)
    ax6.hist(all_steps_corrected, bins=50, alpha=0.5, color='blue', 
             label=f'Drift corrected ({len(all_steps_corrected)} steps)', density=True)
    ax6.set_xlabel('Step size (Å)')
    ax6.set_ylabel('Probability density')
    ax6.set_title(f'Y-direction Step Size Distribution\n{n_molecules} molecules combined')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_multiple_molecule_trajectories.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print statistics for all molecules
    print(f"\nMultiple Molecule Statistics ({n_molecules} molecules):")
    print("-" * 80)
    print(f"{'Mol':<4} {'Y_disp_uncorr':<12} {'Y_disp_corr':<12} {'RMS_uncorr':<10} {'RMS_corr':<10} {'D_single':<12} {'R²':<6}")
    print("-" * 80)
    
    for stats in molecule_stats:
        print(f"{stats['idx']:<4} {stats['y_disp_uncorr']:<12.2f} {stats['y_disp_corr']:<12.2f} "
              f"{stats['rms_uncorr']:<10.2f} {stats['rms_corr']:<10.2f} "
              f"{stats['D_single']:<12.3e} {stats['r_squared']:<6.3f}")
    
    print("-" * 80)
    
    # Summary statistics
    valid_D = [s['D_single'] for s in molecule_stats if not np.isnan(s['D_single'])]
    if valid_D:
        print(f"Average D_single: {np.mean(valid_D):.3e} ± {np.std(valid_D):.3e} Å²/ps")
        print(f"D_single range: {np.min(valid_D):.3e} - {np.max(valid_D):.3e} Å²/ps")
    
    avg_step_uncorr = np.sqrt(np.mean(np.array(all_steps_uncorrected)**2))
    avg_step_corr = np.sqrt(np.mean(np.array(all_steps_corrected)**2))
    print(f"Average step RMS (uncorrected): {avg_step_uncorr:.3f} Å")
    print(f"Average step RMS (drift corrected): {avg_step_corr:.3f} Å")
    print(f"Step size reduction factor: {avg_step_uncorr/avg_step_corr:.1f}")
    
    return molecule_stats

def fit_diffusion_coefficient(time, msd, dimension=1, start_fit=0.2, end_fit=0.8):
    """
    Fit diffusion coefficient from MSD with sensitivity analysis
    
    Parameters:
    -----------
    time : np.ndarray
        Time points
    msd : np.ndarray
        Mean squared displacement
    dimension : int
        Dimensionality for diffusion coefficient calculation
    start_fit : float
        Start of fitting region (fraction of total time)
    end_fit : float
        End of fitting region (fraction of total time)
        
    Returns:
    --------
    tuple : (diffusion_coeff, slope, r_squared, fit_sensitivity)
    """
    start_idx = int(start_fit * len(time))
    end_idx = int(end_fit * len(time))
    
    # Main fit
    slope, intercept, r_value, _, _ = stats.linregress(
        time[start_idx:end_idx], 
        msd[start_idx:end_idx]
    )
    
    diffusion_coeff = slope / (2 * dimension)
    r_squared = r_value**2
    
    # Sensitivity analysis - test different fitting windows
    sensitivity_results = []
    for start_test in [0.1, 0.15, 0.25, 0.3]:
        for end_test in [0.7, 0.75, 0.85, 0.9]:
            if start_test < end_test:
                start_test_idx = int(start_test * len(time))
                end_test_idx = int(end_test * len(time))
                
                if end_test_idx - start_test_idx > 10:  # Minimum points for fit
                    slope_test, _, r_test, _, _ = stats.linregress(
                        time[start_test_idx:end_test_idx],
                        msd[start_test_idx:end_test_idx]
                    )
                    D_test = slope_test / (2 * dimension)
                    sensitivity_results.append(D_test)
    
    fit_sensitivity = np.std(sensitivity_results) / np.mean(sensitivity_results) if sensitivity_results else 0.0
    
    return diffusion_coeff, slope, r_squared, fit_sensitivity

def bootstrap_diffusion_analysis(time, msd, n_bootstrap=1000, start_fit=0.2, end_fit=0.8):
    """
    Bootstrap analysis for uncertainty quantification
    
    Parameters:
    -----------
    time : np.ndarray
        Time points
    msd : np.ndarray
        Mean squared displacement
    n_bootstrap : int
        Number of bootstrap samples
    start_fit : float
        Start of fitting region
    end_fit : float
        End of fitting region
        
    Returns:
    --------
    dict : Bootstrap statistics
    """
    print(f"Performing bootstrap analysis ({n_bootstrap} samples)...")
    
    bootstrap_D = []
    start_idx = int(start_fit * len(time))
    end_idx = int(end_fit * len(time))
    
    time_fit = time[start_idx:end_idx]
    msd_fit = msd[start_idx:end_idx]
    
    for _ in range(n_bootstrap):
        # Resample data points with replacement
        indices = np.random.choice(len(time_fit), size=len(time_fit), replace=True)
        time_boot = time_fit[indices]
        msd_boot = msd_fit[indices]
        
        # Sort by time (important for regression)
        sort_idx = np.argsort(time_boot)
        time_boot = time_boot[sort_idx]
        msd_boot = msd_boot[sort_idx]
        
        # Fit and calculate D
        try:
            slope, _, _, _, _ = stats.linregress(time_boot, msd_boot)
            D_boot = slope / 2  # 1D diffusion
            bootstrap_D.append(D_boot)
        except:
            continue  # Skip failed fits
    
    bootstrap_D = np.array(bootstrap_D)
    
    # Calculate statistics
    mean_D = np.mean(bootstrap_D)
    std_D = np.std(bootstrap_D)
    ci_lower = np.percentile(bootstrap_D, 2.5)
    ci_upper = np.percentile(bootstrap_D, 97.5)
    
    return {
        'mean': mean_D,
        'std': std_D,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'samples': bootstrap_D
    }

def create_comprehensive_plots(time, msd_x, msd_y, msd_z, msd_total, fit_results, 
                             bootstrap_results=None, output_prefix="diffusion_results"):
    """
    Create comprehensive visualization of diffusion analysis
    """
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. MSD components (linear scale)
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(time, msd_x, 'b-', label='MSD_x', alpha=0.7)
    ax1.plot(time, msd_y, 'g-', label='MSD_y', alpha=0.7)
    ax1.plot(time, msd_z, 'r-', label='MSD_z', alpha=0.7, linewidth=2)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('MSD (Å²)')
    ax1.set_title('Mean Squared Displacement Components')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. MSD Z-direction with fit
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(time, msd_z, 'r-', label='MSD_z', alpha=0.7)
    
    D_z, slope, r2, sensitivity = fit_results
    start_idx = int(0.2 * len(time))
    end_idx = int(0.8 * len(time))
    fit_line = slope * time + (msd_z[start_idx] - slope * time[start_idx])
    ax2.plot(time[start_idx:end_idx], fit_line[start_idx:end_idx], 
             'k--', linewidth=2, 
             label=f'Fit: D_z = {D_z:.3e} Å²/ps\\nR² = {r2:.3f}\\nSensitivity = {sensitivity:.1%}')
    
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('MSD_z (Å²)')
    ax2.set_title('Z-Direction Diffusion with Linear Fit')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Log-log plot for diffusive regime validation
    ax3 = plt.subplot(2, 3, 3)
    ax3.loglog(time[1:], msd_z[1:], 'r-', alpha=0.7)
    ax3.loglog(time[1:], time[1:] * msd_z[1] / time[1], 'k--', alpha=0.5, label='Slope = 1')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('MSD_z (Å²)')
    ax3.set_title('Log-Log Plot (Diffusive Regime Check)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Bootstrap distribution (if available)
    if bootstrap_results:
        ax4 = plt.subplot(2, 3, 4)
        ax4.hist(bootstrap_results['samples'], bins=50, alpha=0.7, color='skyblue', density=True)
        ax4.axvline(bootstrap_results['mean'], color='red', linestyle='--', 
                   label=f"Mean: {bootstrap_results['mean']:.3e}")
        ax4.axvline(bootstrap_results['ci_lower'], color='orange', linestyle=':', 
                   label=f"95% CI")
        ax4.axvline(bootstrap_results['ci_upper'], color='orange', linestyle=':')
        ax4.set_xlabel('D_z (Å²/ps)')
        ax4.set_ylabel('Probability Density')
        ax4.set_title('Bootstrap Distribution')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    # 5. Anisotropy analysis
    ax5 = plt.subplot(2, 3, 5)
    anisotropy = msd_z / ((msd_x + msd_y) / 2)
    ax5.plot(time, anisotropy, 'purple', alpha=0.7)
    ax5.axhline(y=1, color='k', linestyle='--', alpha=0.5, label='Isotropic')
    ax5.set_xlabel('Time (ps)')
    ax5.set_ylabel('Anisotropy (D_z/D_xy)')
    ax5.set_title('Diffusion Anisotropy')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Fitting sensitivity analysis
    ax6 = plt.subplot(2, 3, 6)
    # Test different fitting windows
    start_fractions = np.linspace(0.1, 0.3, 10)
    end_fractions = np.linspace(0.7, 0.9, 10)
    D_values = []
    
    for start_frac in start_fractions:
        for end_frac in end_fractions:
            if start_frac < end_frac:
                start_idx = int(start_frac * len(time))
                end_idx = int(end_frac * len(time))
                if end_idx - start_idx > 10:
                    slope_test, _, _, _, _ = stats.linregress(time[start_idx:end_idx], msd_z[start_idx:end_idx])
                    D_values.append(slope_test / 2)
    
    if D_values:
        ax6.hist(D_values, bins=20, alpha=0.7, color='lightcoral')
        ax6.axvline(D_z, color='red', linestyle='--', label=f'Selected fit: {D_z:.3e}')
        ax6.set_xlabel('D_z (Å²/ps)')
        ax6.set_ylabel('Frequency')
        ax6.set_title('Fitting Window Sensitivity')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create separate publication-quality plot for main result
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(time, msd_z, 'r-', label='MSD_z', alpha=0.8, linewidth=1.5)
    
    start_idx = int(0.2 * len(time))
    end_idx = int(0.8 * len(time))
    fit_line = slope * time + (msd_z[start_idx] - slope * time[start_idx])
    ax.plot(time[start_idx:end_idx], fit_line[start_idx:end_idx], 
            'k--', linewidth=2, label='Linear Fit')
    
    # Add bootstrap confidence interval if available
    if bootstrap_results:
        ax.set_title(f'Z-Direction Diffusion Coefficient\n' +
                    f'D_z = {bootstrap_results["mean"]:.3e} ± {bootstrap_results["std"]:.3e} Å²/ps\n' +
                    f'95% CI: [{bootstrap_results["ci_lower"]:.3e}, {bootstrap_results["ci_upper"]:.3e}]')
    else:
        ax.set_title(f'Z-Direction Diffusion Coefficient\nD_z = {D_z:.3e} Å²/ps (R² = {r2:.3f})')
    
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('MSD_z (Å²)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_main_result.png', dpi=300, bbox_inches='tight')
    plt.close()

def apply_finite_size_correction(D_uncorrected, temperature, viscosity, box_length):
    """
    Apply Yeh-Hummer finite size correction
    
    Parameters:
    -----------
    D_uncorrected : float
        Uncorrected diffusion coefficient (Å²/ps)
    temperature : float
        Temperature (K)
    viscosity : float
        Viscosity (Pa·s)
    box_length : float
        Box length (Å)
        
    Returns:
    --------
    tuple : (D_corrected, correction)
    """
    # Convert units for calculation
    box_length_m = box_length * 1e-10  # Å to m
    D_uncorrected_m2s = D_uncorrected * 1e-8  # Å²/ps to m²/s
    
    # Yeh-Hummer correction
    correction_m2s = XI * KB * temperature / (6 * np.pi * viscosity * box_length_m)
    
    D_corrected_m2s = D_uncorrected_m2s + correction_m2s
    
    # Convert back to Å²/ps
    correction = correction_m2s * 1e8
    D_corrected = D_corrected_m2s * 1e8
    
    return D_corrected, correction

def main():
    parser = argparse.ArgumentParser(
        description='Calculate diffusion coefficients following best practices',
        epilog='Note: Input trajectory MUST be unwrapped using: gmx trjconv -mol -nojump'
    )
    parser.add_argument('topology', help='Topology file (.gro, .pdb)')
    parser.add_argument('trajectory', help='UNWRAPPED trajectory file (.xtc, .trr)')
    parser.add_argument('--selection', default='all', help='Atom selection (default: all)')
    parser.add_argument('--fragment-size', type=int, default=28, 
                       help='Atoms per molecular fragment (default: 28)')
    parser.add_argument('--temperature', type=float, default=300, 
                       help='System temperature (default: 300 K)')
    parser.add_argument('--viscosity', type=float, 
                       help='System viscosity (Pa·s) for finite size correction')
    parser.add_argument('--bootstrap', action='store_true', 
                       help='Perform bootstrap uncertainty analysis')
    parser.add_argument('--n-bootstrap', type=int, default=1000, 
                       help='Number of bootstrap samples (default: 1000)')
    parser.add_argument('--start-fit', type=float, default=0.2, 
                       help='Start of fitting region (fraction, default: 0.2)')
    parser.add_argument('--end-fit', type=float, default=0.8, 
                       help='End of fitting region (fraction, default: 0.8)')
    parser.add_argument('--output', default='diffusion_results', 
                       help='Output prefix (default: diffusion_results)')
    parser.add_argument('--validate-drift', action='store_true',
                       help='Create drift correction validation plots')
    parser.add_argument('--molecule-indices', type=int, nargs='+', default=[0],
                       help='Indices of molecules for trajectory analysis (default: [0])')
    
    args = parser.parse_args()
    
    print("="*60)
    print("DIFFUSION COEFFICIENT CALCULATION")
    print("Following LiveCoMS Best Practices (Maginn et al., 2019)")
    print("="*60)
    
    # Load universe
    print(f"\nLoading trajectory: {args.trajectory}")
    print(f"Using topology: {args.topology}")
    u = mda.Universe(args.topology, args.trajectory)
    
    # Validate inputs
    if not validate_unwrapped_trajectory(u, args.fragment_size):
        print("\nERROR: Trajectory validation failed!")
        print("Please unwrap trajectory using: gmx trjconv -mol -nojump")
        return 1
    
    validate_simulation_quality(u)
    
    # Compute system drift
    drift = compute_system_drift(u)
    
    # Compute fragment deviations
    print(f"\nAnalyzing molecular fragments...")
    deviations, fragments = compute_fragment_deviations(u, args.selection, args.fragment_size)
    
    # Validation plots if requested
    if args.validate_drift:
        plot_drift_correction_validation(
            np.array([ts.time for ts in u.trajectory]), 
            deviations, drift, args.output
        )
    
    # Multiple molecule trajectories if requested
    if args.molecule_indices:
        # Filter valid molecule indices
        valid_indices = [idx for idx in args.molecule_indices if idx < deviations.shape[1]]
        if len(valid_indices) != len(args.molecule_indices):
            invalid_indices = [idx for idx in args.molecule_indices if idx >= deviations.shape[1]]
            print(f"Warning: Invalid molecule indices {invalid_indices} (max: {deviations.shape[1]-1})")
        
        if valid_indices:
            molecule_stats = plot_multiple_molecule_trajectories(
                deviations, drift, 
                np.array([ts.time for ts in u.trajectory]), 
                args.output, valid_indices
            )
    
    # Calculate MSD components - NOTE: Intentionally swapped y and z for y-axis analysis
    print("Calculating mean squared displacements...")
    msd_x, msd_z, msd_y, msd_total = compute_msd_components(deviations, drift)
    
    # Get time array
    time = np.array([ts.time for ts in u.trajectory])
    
    # Validate diffusive regime - using msd_z which contains y-axis data
    box_length = u.dimensions[1]  # Y-dimension for cylinder axis
    if not validate_diffusive_regime(msd_z, time, box_length):
        warnings.warn("Diffusive regime validation failed - results may be unreliable")
    
    # Fit diffusion coefficient
    print(f"\nFitting diffusion coefficient (fitting region: {args.start_fit:.1%} - {args.end_fit:.1%})...")
    D_z, slope, r2, sensitivity = fit_diffusion_coefficient(
        time, msd_z, dimension=1, start_fit=args.start_fit, end_fit=args.end_fit
    )
    
    # Bootstrap analysis if requested
    bootstrap_results = None
    if args.bootstrap:
        bootstrap_results = bootstrap_diffusion_analysis(
            time, msd_z, args.n_bootstrap, args.start_fit, args.end_fit
        )
    
    # Apply finite size correction if viscosity provided
    D_z_corrected = None
    correction = None
    if args.viscosity:
        print(f"\nApplying finite size correction...")
        print(f"Using viscosity: {args.viscosity:.6f} Pa·s")
        D_z_corrected, correction = apply_finite_size_correction(
            D_z, args.temperature, args.viscosity, box_length
        )
    
    # Create comprehensive plots
    print(f"\nGenerating plots...")
    create_comprehensive_plots(
        time, msd_x, msd_y, msd_z, msd_total, 
        (D_z, slope, r2, sensitivity), bootstrap_results, args.output
    )
    
    # Save numerical results
    print(f"Saving numerical data...")
    np.savetxt(f'{args.output}_msd_data.txt', 
               np.column_stack((time, msd_x, msd_y, msd_z, msd_total)),
               header='Time(ps) MSD_x(Å²) MSD_y(Å²) MSD_z(Å²) MSD_total(Å²)',
               fmt='%.6f')
    
    np.savetxt(f'{args.output}_drift.txt', drift,
               header='Drift_x(Å) Drift_y(Å) Drift_z(Å)',
               fmt='%.6f')
    
    # Print comprehensive results
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    
    print(f"\nSimulation Parameters:")
    print(f"  Temperature: {args.temperature:.1f} K")
    print(f"  Box length (Y-axis): {box_length:.2f} Å")
    print(f"  Number of fragments: {len(fragments)}")
    print(f"  Trajectory length: {len(time)} frames ({time[-1]:.1f} ps)")
    
    print(f"\nDiffusion Coefficient (Y-direction, stored as msd_z):")
    print(f"  D_y = {D_z:.6e} Å²/ps")
    print(f"  D_y = {D_z*1e4:.6f} nm²/μs")
    print(f"  Linear fit R² = {r2:.4f}")
    print(f"  Fitting sensitivity = {sensitivity:.1%}")
    
    if bootstrap_results:
        print(f"\nBootstrap Uncertainty Analysis:")
        print(f"  Mean D_y = {bootstrap_results['mean']:.6e} ± {bootstrap_results['std']:.6e} Å²/ps")
        print(f"  95% Confidence Interval: [{bootstrap_results['ci_lower']:.6e}, {bootstrap_results['ci_upper']:.6e}] Å²/ps")
        print(f"  Relative uncertainty = {bootstrap_results['std']/bootstrap_results['mean']:.1%}")
    
    if D_z_corrected is not None:
        print(f"\nFinite Size Correction (Yeh-Hummer):")
        print(f"  Correction = +{correction:.6e} Å²/ps")
        print(f"  D_y (corrected) = {D_z_corrected:.6e} Å²/ps")
        print(f"  D_y (corrected) = {D_z_corrected*1e4:.6f} nm²/μs")
        print(f"  Relative correction = {correction/D_z:.1%}")
    
    print(f"\nAnisotropy Analysis:")
    final_D_xz = (msd_x[-1] + msd_y[-1]) / (4 * time[-1])  # Average xz diffusion
    anisotropy = D_z / final_D_xz if final_D_xz > 0 else float('inf')
    print(f"  D_xz (estimated) = {final_D_xz:.6e} Å²/ps")
    print(f"  Anisotropy (D_y/D_xz) = {anisotropy:.2f}")
    
    print(f"\nQuality Metrics:")
    print(f"  System drift magnitude = {np.mean(np.linalg.norm(drift, axis=1)):.4f} Å/frame")
    print(f"  Final RMS displacement = {np.sqrt(msd_z[-1]):.2f} Å")
    print(f"  RMS/box ratio = {np.sqrt(msd_z[-1])/box_length:.3f}")
    
    # Save comprehensive results to file
    with open(f'{args.output}_summary.txt', 'w') as f:
        f.write("Diffusion Coefficient Analysis Summary\n")
        f.write("="*50 + "\n\n")
        
        f.write(f"Input Files:\n")
        f.write(f"  Topology: {args.topology}\n")
        f.write(f"  Trajectory: {args.trajectory}\n")
        f.write(f"  Selection: {args.selection}\n\n")
        
        f.write(f"Analysis Parameters:\n")
        f.write(f"  Fragment size: {args.fragment_size} atoms\n")
        f.write(f"  Temperature: {args.temperature} K\n")
        f.write(f"  Fitting region: {args.start_fit:.1%} - {args.end_fit:.1%}\n")
        if args.viscosity:
            f.write(f"  System viscosity: {args.viscosity:.6f} Pa·s\n")
        f.write(f"\n")
        
        f.write(f"System Properties:\n")
        f.write(f"  Box length (Y-axis): {box_length:.2f} Å")
        f.write(f"  Box dimensions: {u.dimensions[0]:.1f} x {u.dimensions[1]:.1f} x {u.dimensions[2]:.1f} Å\n")
        f.write(f"  Number of fragments: {len(fragments)}\n")
        f.write(f"  Trajectory length: {len(time)} frames ({time[-1]:.1f} ps)\n\n")
        
        f.write(f"Diffusion Coefficient Results (Y-direction):\n")
        f.write(f"  D_y = {D_z:.6e} Å²/ps = {D_z*1e4:.6f} nm²/μs\n")
        f.write(f"  Linear fit R² = {r2:.4f}\n")
        f.write(f"  Fitting sensitivity = {sensitivity:.1%}\n")
        
        if bootstrap_results:
            f.write(f"\n  Bootstrap Analysis:\n")
            f.write(f"    Mean ± Std = {bootstrap_results['mean']:.6e} ± {bootstrap_results['std']:.6e} Å²/ps\n")
            f.write(f"    95% CI = [{bootstrap_results['ci_lower']:.6e}, {bootstrap_results['ci_upper']:.6e}] Å²/ps\n")
            f.write(f"    Relative uncertainty = {bootstrap_results['std']/bootstrap_results['mean']:.1%}\n")
        
        if D_z_corrected is not None:
            f.write(f"\n  Finite Size Correction:\n")
            f.write(f"    Correction = +{correction:.6e} Å²/ps ({correction/D_z:.1%})\n")
            f.write(f"    D_y (corrected) = {D_z_corrected:.6e} Å²/ps = {D_z_corrected*1e4:.6f} nm²/μs\n")
        
        f.write(f"\nAnisotropy:\n")
        f.write(f"  D_xz (estimated) = {final_D_xz:.6e} Å²/ps\n")
        f.write(f"  Anisotropy (D_y/D_xz) = {anisotropy:.2f}\n")
        
        f.write(f"\nQuality Metrics:\n")
        f.write(f"  System drift = {np.mean(np.linalg.norm(drift, axis=1)):.4f} Å/frame\n")
        f.write(f"  Final RMS displacement = {np.sqrt(msd_z[-1]):.2f} Å\n")
        f.write(f"  RMS/box ratio = {np.sqrt(msd_z[-1])/box_length:.3f}\n")
    
    # Recommendations and warnings
    print(f"\n" + "="*60)
    print("RECOMMENDATIONS AND WARNINGS")
    print("="*60)
    
    if r2 < 0.95:
        print(f"⚠️  Low R² ({r2:.3f}) suggests poor linear fit - consider different fitting region")
    
    if sensitivity > 0.1:
        print(f"⚠️  High fitting sensitivity ({sensitivity:.1%}) - results may depend on fitting window choice")
    
    if np.sqrt(msd_z[-1])/box_length < 0.2:
        print(f"⚠️  Small displacement relative to box size - consider longer simulation")
    
    # Special warnings for high drift systems
    drift_magnitude = np.mean(np.linalg.norm(drift, axis=1))
    if drift_magnitude > 1.0:
        print(f"⚠️  Large system drift ({drift_magnitude:.3f} Å/frame) detected")
        print(f"    Consider using --validate-drift to check drift correction effectiveness")
    
    if not args.viscosity:
        print(f"ℹ️  Finite size correction not applied - provide --viscosity for more accurate results")
    
    if not args.bootstrap:
        print(f"ℹ️  Bootstrap analysis not performed - use --bootstrap for uncertainty quantification")
    
    if anisotropy < 0.1:
        print(f"ℹ️  Very low anisotropy (D_y/D_xz = {anisotropy:.3f}) - may indicate confined diffusion")
    elif anisotropy > 10:
        print(f"ℹ️  High anisotropy (D_y/D_xz = {anisotropy:.1f}) - fast axial diffusion confirmed")
    
    # Compare to expected values
    expected_D = 0.8e-4  # nm²/μs to Å²/ps conversion: 0.8 * 1e-4
    if D_z > 0:
        ratio_to_expected = (D_z * 1e4) / 0.8  # Convert to nm²/μs and compare
        print(f"ℹ️  D_y = {D_z*1e4:.4f} nm²/μs vs expected ~0.8 nm²/μs (ratio: {ratio_to_expected:.2f})")
    
    print(f"\n✅ Analysis complete! Files saved with prefix: {args.output}")
    print(f"📊 Main result plot: {args.output}_main_result.png")
    print(f"📈 Comprehensive analysis: {args.output}_comprehensive_analysis.png")
    if args.validate_drift:
        print(f"🔍 Drift validation: {args.output}_drift_validation.png")
    if args.molecule_indices:
        print(f"🎯 Multiple molecules: {args.output}_multiple_molecule_trajectories.png")
    print(f"📄 Summary report: {args.output}_summary.txt")
    
    return 0

if __name__ == '__main__':
    import sys
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)