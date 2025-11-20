#!/bin/bash

# Create necessary directories if they don't exist
mkdir -p "../figures/"
mkdir -p "../tmp"

echo "Processing equilibration data from eq1 to eq5..."

# Define the output filename
OUTPUT_NAME="eq_plot"

# Create a temporary directory to store extracted data
mkdir -p "../tmp/${OUTPUT_NAME}"

# Loop through eq1 to eq5 and extract energy data
for i in {1..5}; do
    echo "Extracting data from eq${i}..."
    
    # Extract temperature, z-coordinate, and density data
    printf '17\n0\n' | gmx energy -f "../CG-reference/eq${i}/eq.edr" -o "../tmp/${OUTPUT_NAME}/eq${i}_z.xvg"
    printf '19\n0\n' | gmx energy -f "../CG-reference/eq${i}/eq.edr" -o "../tmp/${OUTPUT_NAME}/eq${i}_density.xvg"
done

# Extract data from the production run
echo "Extracting data from production run..."
printf '17\n0\n' | gmx energy -f "../CG-reference/run/run.edr" -o "../tmp/${OUTPUT_NAME}/run_z.xvg"
printf '19\n0\n' | gmx energy -f "../CG-reference/run/run.edr" -o "../tmp/${OUTPUT_NAME}/run_density.xvg"

# Create Python script to plot the data
cat > plot_eq_script.py << 'EOL'
import matplotlib.pyplot as plt
import numpy as np
import os

# Function to read XVG files (GROMACS format)
def read_xvg(filename):
    time = []
    value = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('@'):
                    continue
                cols = line.strip().split()
                if len(cols) >= 2:
                    time.append(float(cols[0]))
                    value.append(float(cols[1]))
    except FileNotFoundError:
        print(f"Warning: File {filename} not found, returning empty arrays")
    return np.array(time), np.array(value)

# Function to concatenate multiple datasets with proper time offset
def concatenate_datasets(time_values_pairs):
    if not time_values_pairs:
        return np.array([]), np.array([])
    
    # Initialize with the first dataset
    concatenated_time, concatenated_values = time_values_pairs[0]
    
    # Add subsequent datasets with proper time offsets
    for i in range(1, len(time_values_pairs)):
        current_time, current_values = time_values_pairs[i]
        
        if len(current_time) == 0:
            continue
            
        if len(concatenated_time) == 0:
            concatenated_time = current_time
            concatenated_values = current_values
            continue
        
        # Calculate time offset (next dataset should start where previous ends)
        time_offset = concatenated_time[-1]
        
        # Adjust time by adding the offset
        adjusted_time = current_time + time_offset
        
        # Concatenate the arrays
        concatenated_time = np.concatenate([concatenated_time, adjusted_time])
        concatenated_values = np.concatenate([concatenated_values, current_values])
    
    return concatenated_time, concatenated_values

try:
    # Read data for eq1 through eq5 and production run
    datasets = {}
    transition_times = []
    transition_labels = []
    last_time = 0
    
    for metric in ["z", "density"]:
        time_values_pairs = []
        
        # Read equilibration data
        for i in range(1, 6):  # eq1 to eq5
            filename = f"../tmp/eq_plot/eq{i}_{metric}.xvg"
            time, values = read_xvg(filename)
            
            # Only add non-empty datasets
            if len(time) > 0:
                time_values_pairs.append((time, values))
                
                # Store transition time
                if metric == "temp" and i > 1:
                    transition_times.append(last_time)
                    transition_labels.append(f"eq{i-1} → eq{i}")
                
                if metric == "temp" and len(time) > 0:
                    last_time += time[-1]
        
        # Add production run data
        filename = f"../tmp/eq_plot/run_{metric}.xvg"
        time, values = read_xvg(filename)
        
        if len(time) > 0:
            time_values_pairs.append((time, values))
            
            # Store transition time between eq5 and production run
            if metric == "temp":
                transition_times.append(last_time)
                transition_labels.append("eq5 → production")
        
        # Concatenate all datasets for this metric
        datasets[metric] = concatenate_datasets(time_values_pairs)
    
    # Create a plot showing the continuous combined data
    fig, ax1 = plt.subplots(figsize=(14, 8))
    
    # Temperature plot (primary y-axis)
    color = 'tab:red'
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Density (kg/m³)', color=color)
    line1 = ax1.plot(datasets["density"][0], datasets["density"][1], color=color, label='Density')
    ax1.tick_params(axis='y', labelcolor=color)
    
    # Z-coordinate plot (secondary y-axis)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Z-coordinate (nm)', color=color)
    line2 = ax2.plot(datasets["z"][0], datasets["z"][1], color=color, label='Z-coordinate')
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Add vertical lines to mark transitions between steps
    vlines = []
    for i, time in enumerate(transition_times):
        vline = plt.axvline(x=time, color='black', linestyle=':', alpha=0.7)
        if i == 0:  # Only add one line to the legend
            vline.set_label('Step transition')
            vlines.append(vline)

    # Add legend
    lines = line1 + line2  + vlines
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='best')
    
    plt.title('MD Simulation Analysis: Continuous Plot (eq1-eq5 + 10 μs production run)')
    plt.tight_layout()
    
    # Add annotations for the transitions
    y_pos = ax1.get_ylim()[1] * 0.95  # Position annotations near the top
    for i, (time, label) in enumerate(zip(transition_times, transition_labels)):
        ax1.annotate(label, (time, y_pos), 
                   xytext=(0, 10), textcoords='offset points',
                   rotation=90, ha='right', va='top', fontsize=8)
    
    # Save figure
    plt.savefig("../figures/eq_plot.png", dpi=300)
    print("Plot saved to ../figures/eq_plot.png")
    
except Exception as e:
    print(f"Error plotting data: {e}")
    import traceback
    traceback.print_exc()
EOL

# Run the Python script
python3 plot_eq_script.py

# Clean up
echo "Cleaning up temporary files..."
rm plot_eq_script.py

echo "Done!"