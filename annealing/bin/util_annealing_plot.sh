#!/bin/bash

# Check if OUTPUT_NAME is provided
if [ $# -eq 0 ]; then
    echo "Please provide OUTPUT_NAME as an argument"
    echo "Usage: $0 <OUTPUT_NAME>"
    exit 1
fi

mkdir -p ../figures/

OUTPUT_NAME=$1
echo "Processing for OUTPUT_NAME: $OUTPUT_NAME"

# Create necessary directories if they don't exist
mkdir -p "../system/${OUTPUT_NAME}"
mkdir -p "../tmp"

# Run the gmx energy commands
printf '15\n0\n' | gmx energy -f "../system/${OUTPUT_NAME}/annealing_${OUTPUT_NAME}.edr" -o "../system/${OUTPUT_NAME}/annealing_${OUTPUT_NAME}_temp.xvg"
printf '21\n0\n' | gmx energy -f "../system/${OUTPUT_NAME}/annealing_${OUTPUT_NAME}.edr" -o "../system/${OUTPUT_NAME}/annealing_${OUTPUT_NAME}_z.xvg"
printf '23\n0\n' | gmx energy -f "../system/${OUTPUT_NAME}/annealing_${OUTPUT_NAME}.edr" -o "../system/${OUTPUT_NAME}/annealing_${OUTPUT_NAME}_density.xvg"

echo "Extracting NPT simulation data..."
printf '15\n0\n' | gmx energy -f "../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}.edr" -o "../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}_temp.xvg"
printf '21\n0\n' | gmx energy -f "../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}.edr" -o "../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}_z.xvg"
printf '23\n0\n' | gmx energy -f "../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}.edr" -o "../system/${OUTPUT_NAME}/annealing_npt_${OUTPUT_NAME}_density.xvg"

# Now run the Python script to plot the data
cat > plot_script.py << 'EOL'
import sys
import matplotlib.pyplot as plt
import numpy as np
import os

# Get OUTPUT_NAME from command line argument
output_name = sys.argv[1]

# Define file paths for annealing simulation
annealing_temp_file = f"../system/{output_name}/annealing_{output_name}_temp.xvg"
annealing_z_file = f"../system/{output_name}/annealing_{output_name}_z.xvg"
annealing_density_file = f"../system/{output_name}/annealing_{output_name}_density.xvg"

# Define file paths for NPT simulation
npt_temp_file = f"../system/{output_name}/annealing_npt_{output_name}_temp.xvg"
npt_z_file = f"../system/{output_name}/annealing_npt_{output_name}_z.xvg"
npt_density_file = f"../system/{output_name}/annealing_npt_{output_name}_density.xvg"

# Function to read XVG files (GROMACS format)
def read_xvg(filename):
    time = []
    value = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            cols = line.strip().split()
            if len(cols) >= 2:
                time.append(float(cols[0]))
                value.append(float(cols[1]))
    return np.array(time), np.array(value)

# Function to concatenate two datasets with proper time offset
def concatenate_data(time1, values1, time2, values2):
    # If first dataset is empty, just return the second dataset
    if len(time1) == 0:
        return time2, values2
    
    # If second dataset is empty, just return the first dataset
    if len(time2) == 0:
        return time1, values1
    
    # Calculate time offset (time2 should start where time1 ends)
    time_offset = time1[-1]
    
    # Adjust time2 by adding the offset
    adjusted_time2 = time2 + time_offset
    
    # Concatenate the arrays
    concatenated_time = np.concatenate([time1, adjusted_time2])
    concatenated_values = np.concatenate([values1, values2])
    
    return concatenated_time, concatenated_values

try:
    # Read annealing data
    annealing_temp_time, annealing_temp_values = read_xvg(annealing_temp_file)
    annealing_z_time, annealing_z_values = read_xvg(annealing_z_file)
    annealing_density_time, annealing_density_values = read_xvg(annealing_density_file)
    
    # Read NPT data
    npt_temp_time, npt_temp_values = read_xvg(npt_temp_file)
    npt_z_time, npt_z_values = read_xvg(npt_z_file)
    npt_density_time, npt_density_values = read_xvg(npt_density_file)
    
    # Concatenate annealing and NPT data
    combined_temp_time, combined_temp_values = concatenate_data(
        annealing_temp_time, annealing_temp_values, npt_temp_time, npt_temp_values)
    
    combined_z_time, combined_z_values = concatenate_data(
        annealing_z_time, annealing_z_values, npt_z_time, npt_z_values)
    
    combined_density_time, combined_density_values = concatenate_data(
        annealing_density_time, annealing_density_values, npt_density_time, npt_density_values)
        
    # Create a second plot showing the continuous combined data
    fig2, ax1 = plt.subplots(figsize=(14, 8))
    
    # Temperature plot (primary y-axis)
    color = 'tab:red'
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)', color=color)
    line1 = ax1.plot(combined_temp_time, combined_temp_values, color=color, label='Temperature')
    ax1.tick_params(axis='y', labelcolor=color)
    
    # Z-coordinate plot (secondary y-axis)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Z-coordinate (nm)', color=color)
    line2 = ax2.plot(combined_z_time, combined_z_values, color=color, label='Z-coordinate')
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Density plot (third y-axis)
    ax3 = ax1.twinx()
    # Offset the 3rd axis
    ax3.spines["right"].set_position(("axes", 1.1))
    color = 'tab:green'
    ax3.set_ylabel('Density (kg/m³)', color=color)
    line3 = ax3.plot(combined_density_time, combined_density_values, color=color, label='Density')
    ax3.tick_params(axis='y', labelcolor=color)
    
    # Add vertical line to mark the transition between annealing and NPT
    if len(annealing_temp_time) > 0:
        transition_time = annealing_temp_time[-1]
        plt.axvline(x=transition_time, color='black', linestyle=':', 
                   label='Annealing → NPT transition')
    
    # Add legend
    lines = line1 + line2 + line3 + [plt.Line2D([0], [0], color='black', linestyle=':')]
    labels = [l.get_label() for l in lines[:-1]] + ['Annealing → NPT transition']
    ax1.legend(lines, labels, loc='best')
    
    plt.title(f'MD Simulation Analysis: Continuous Plot for {output_name}')
    plt.tight_layout()
    
    # Save figure
    plt.savefig(f"../figures/{output_name}_continuous_annealing_npt.png", dpi=300)
    print(f"Plot saved to ../figures/{output_name}_continuous_annealing_npt.png")
    
except Exception as e:
    print(f"Error plotting data: {e}")
    import traceback
    traceback.print_exc()
EOL

# Run the Python script with the OUTPUT_NAME as an argument
python3 plot_script.py "$OUTPUT_NAME"

# Remove the temporary script
rm plot_script.py
