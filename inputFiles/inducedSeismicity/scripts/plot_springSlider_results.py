import argparse
import numpy as np
from geos.hdf5_wrapper import hdf5_wrapper
import matplotlib.pyplot as plt
import os

import numpy as np
import h5py  # Assuming hdf5_wrapper is a wrapper around h5py

def remove_padding(data):
    """Removes trailing zeros along the first axis while keeping initial zeros."""
    if data.ndim < 1:
        raise ValueError("remove_padding requires an array with at least one dimension.")

    # Find the last index along axis 0 that contains nonzero elements
    nonzero_mask = np.any(data != 0, axis=tuple(range(1, data.ndim)))
    nonzero_indices = np.where(nonzero_mask)[0]

    if nonzero_indices.size == 0:  # If all elements are zero
        return np.empty((0, *data.shape[1:]), dtype=data.dtype)

    last_nonzero = nonzero_indices[-1]
    return data[:last_nonzero + 1]  # Keep everything up to the last nonzero entry

def getDataFromHDF5(hdf5FilePath, var_name):
    """Extracts data and time from an HDF5 file and removes trailing padding from the data."""
    with h5py.File(hdf5FilePath, 'r') as data:
        var = np.asarray(data[f'{var_name}'])
        time = np.asarray(data[f'{var_name} Time'])

    # Remove trailing zeros from var
    var = remove_padding(var)
    time = time[:var.shape[0]]  # Ensure time matches the new length of var

    return time, var    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, help='Path to hdf5 file')
    args = parser.parse_args()

    normalized_dir = os.path.abspath(args.dir) 
    if not os.path.isdir(normalized_dir):
        print(f"Error: {normalized_dir} is not a valid directory.")
        exit(1)
    
    filePath = f"{normalized_dir}/springSlider.hdf5"
    if not os.path.exists(filePath):
        print(f"Error: {filePath} not found.")
        exit(1)    
    
    time, slipRate = getDataFromHDF5( filePath, "slipRate" )
    time, stateVariable = getDataFromHDF5( filePath, "stateVariable" )
    time, slip = getDataFromHDF5( filePath, "totalSlip" )
    slip_x, slip_y = slip.squeeze(1).T 
    slip_norm = np.sqrt(slip_x**2 + slip_y**2)
    time, deltaSlip = getDataFromHDF5( filePath, "deltaSlip" )
    deltaSlip_x, deltaSlip_y = deltaSlip.squeeze(1).T
    deltaSlip_norm = np.sqrt(deltaSlip_x**2 + deltaSlip_y**2)

    # Convert time to years
    time_in_years = time / (365 * 24 * 3600)  # Assuming time is in seconds

    # Plotting
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8), sharex=True)

    # Plot pressure 
    ax1.set_xlabel('Time (years)')
    ax1.set_ylabel('slip rate (Pa)', color='tab:blue')
    ax1.plot(time_in_years, slipRate, label="slip rate", color='tab:blue')
    ax1.set_yscale('log')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # Set x-axis limits to 0 to 2 years
    ax1.set_xlim(0, np.max(time_in_years))
  
    # Plot stateVariable 
    ax2.set_ylabel('state variable', color='tab:red')
    ax2.plot(time_in_years, stateVariable, label="state variable", color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    ax2.set_xlim(0, np.max(time_in_years))

    ax3.set_ylabel('slip [m]', color='tab:red')
    ax3.plot(time_in_years, slip_norm, label="slip [m]", color='tab:red')
    ax3.tick_params(axis='y', labelcolor='tab:red')
    ax3.set_xlim(0, np.max(time_in_years))

    ax4.set_ylabel('deltaSlip [m]', color='tab:red')
    ax4.plot(time_in_years, deltaSlip_norm, label="slip [m]", color='tab:red')
    ax4.tick_params(axis='y', labelcolor='tab:red')
    ax4.set_xlim(0, np.max(time_in_years))

    # Add grid and title
    fig.suptitle('Spring slider solution', fontsize=16)
    plt.tight_layout()
    plt.grid()

    # Show plot
    plt.show()