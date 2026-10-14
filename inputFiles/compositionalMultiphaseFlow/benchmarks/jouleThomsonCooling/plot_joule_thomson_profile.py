import h5py
import numpy
import json
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.cm as cm
import typing

def get_selected_times() -> dict:
    """
    Defines the specific simulation times for which data will be extracted and plotted.

    Returns:
        dict: A dictionary mapping time in seconds (int) to a human-readable label (str).
    """
    DAY_2_SEC = 86400
    YEAR_2_SEC = int(DAY_2_SEC * 365.25 + 0.5)
    return {
        DAY_2_SEC: '1 day',
        10*DAY_2_SEC: '10 days',
        20*DAY_2_SEC: '20 days',
        30*DAY_2_SEC: '30 days',
        90*DAY_2_SEC: '90 days',
        YEAR_2_SEC: '1 year',
        15*YEAR_2_SEC: '15 years',
    }

def get_time_indices(time_data: numpy.ndarray, selected_times: dict) -> dict:
    """
    Finds the array indices closest to the requested target times.

    Args:
        time_data (numpy.ndarray): The 1D array of time steps from the simulation data.
        selected_times (dict): The dictionary of target times (keys are in seconds).

    Returns:
        dict: A dictionary mapping the target time (in seconds) to its corresponding index in the data arrays.
    """
    indices = {}
    for time in selected_times:
        # Find the index of the value in time_data that is closest to 'time'
        indices[time] = numpy.absolute(time_data - time).argmin()
    return indices

def read_data(file_name: str, selected_times: dict) -> dict:
    """
    Reads simulation results from an HDF5 file for specified time steps.

    Args:
        file_name (str): The path to the HDF5 file containing the simulation data.
        selected_times (dict): The dictionary defining the target times to extract.

    Returns:
        dict: A nested dictionary containing the extracted data. 
              The top-level keys are the target times in seconds.
              Each value is a dictionary containing arrays for radius ('r'), 
              pressure ('p'), temperature ('t'), and CO2 mass fraction ('y').
    """
    data = {}
    with h5py.File(file_name, 'r') as h5:
        # Assume that all the time arrays are the same, use elementCenter Time as reference
        element_center_time = h5['elementCenter Time'][()].flatten()
        time_indices = get_time_indices(element_center_time, selected_times)

        # Read the full arrays for element centers and properties
        element_center = h5['elementCenter'][()]
        pressure = h5['pressure'][()]
        temperature = h5['temperature'][()]
        phase_comp_fraction = h5['fluid_phaseCompFraction'][()]
        
        # Extract data for the closest matching time step for each selected time
        for time in selected_times:
            index = time_indices[time]
            # Assuming radial geometry, take the first component (x-coordinate) as radius
            r = element_center[index, :, 0]
            # Convert Temperature from Kelvin to Celsius
            t = temperature[index, :] - 273.15
            # Convert Pressure from Pascals to Bar
            p = pressure[index, :] * 1.0e-5
            # Extract CO2 Mass Fraction (assuming it's the 4th component, index 3)
            y = phase_comp_fraction[index, :, 3]
            data[time] = {"r": r, "p": p, "t": t, "y": y}

    return data

def plot_profile_data(item: str, data: dict, selected_times: dict, reference: dict = None):
    """
    Plots the spatial profile of a specified property at different time steps.

    Args:
        item (str): The property to plot ('p' for pressure, 't' for temperature, 'y' for mass fraction).
        data (dict): The dictionary containing the simulation data extracted by read_data().
        selected_times (dict): The dictionary mapping times to labels, used for the legend.
        reference (dict, optional): Reference data to plot alongside simulation results. Defaults to None.
    """
    # Set global font size for matplotlib plots
    plt.rcParams.update({'font.size': 16})

    # Setup Colormap to distinguish different time steps
    cmap = plt.colormaps['jet']
    num_times = len(data.keys())
    colors = [cmap((i+0.5)/num_times) for i in range(num_times)]
    
    # Configuration for different plot types
    figs = {
        "p" : { "t": 'Pressure', "u": 'bar', "l": (45, 75) },
        "t" : { "t": 'Temperature', "u": '°C', "l": (30, 50) },
        "y" : { "t": 'CO2 Mass Fraction', "u": '-' }
    }
    
    # Check if the requested item is valid
    if item not in figs:
        print(f"Warning: Property '{item}' is not configured for plotting. Skipping.")
        return

    fig, ax = plt.subplots(figsize=(9, 5))
    sorted_times = [t for t in sorted(data.keys())]
    
    # Plot data for each time step
    for ti, time in enumerate(sorted_times):
        r = data[time]["r"]
        # Plot simulation data as a solid line
        ax.plot(r, data[time][item], color=colors[ti], linewidth=2)
        
        # Plot reference data as points if available
        if reference and str(time) in reference:
            r0, v0 = reference[str(time)][item]
            ax.plot(r0, v0, color=colors[ti], linewidth=0, linestyle="None", marker="o")

    # Create custom legend handles using Patch to match colors
    custom_handles = []
    for ti, time in enumerate(sorted_times):
        time_label = selected_times[time]
        custom_handles.append(Patch(facecolor=colors[ti], edgecolor='black', label=time_label))
  
    # Configure axes labels, title, and scales
    ax.set_xlabel('r (m)')
    ax.set_ylabel("%s (%s)" % (figs[item]["t"], figs[item]["u"]))
    ax.set_title("%s Profile" % (figs[item]["t"]))
    ax.set_xscale('log')
    ax.set_xlim(0.1, 1000)
    
    # Apply y-axis limits if specified in the configuration
    if "l" in figs[item]:
        ax.set_ylim(figs[item]["l"])
        
    ax.grid(True)
    ax.legend(handles=custom_handles, loc='best', fontsize=10)
    fig.tight_layout()

    plt.show(block=True)

def plot_property(property: str, reference: bool = False):
    # Set the directory containing the data
    geos_directory = '../../../..'
    data_directory = os.path.join(geos_directory, "inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling")
    data_directory = os.path.normpath(data_directory)

    # Get target times and extract simulation data
    selected_times = get_selected_times()
    file_name = os.path.join(data_directory, "jouleThomsonCooling.hdf5")
    
    try:
        data = read_data(file_name, selected_times)
    except FileNotFoundError:
        print(f"Error: Could not find data file '{file_name}'. Please ensure the file exists.")
        exit(1)

    # Load reference data if requested
    reference_data = None
    if reference:
        ref_file_name = os.path.join(data_directory, "oldenburg.json")
        try:
            with open(ref_file_name, mode="r") as f:
                reference_data = json.load(f)
        except FileNotFoundError:
             print(f"Warning: Reference file '{ref_file_name}' not found. Plotting without reference data.")

    # Iterate through requested properties and generate plots
    props_to_plot = property.lower()
    for prop in ['p', 't', 'y']:
        if prop in props_to_plot:
            plot_profile_data(prop, data, selected_times, reference=reference_data)

def plot_pressure():
    plot_property("p")

def plot_temperature():
    plot_property("t", True)

def plot_co2_fraction():
    plot_property("y")

if __name__ == "__main__":
    # Setup argparse for command line arguments
    parser = argparse.ArgumentParser(description="Plot Joule-Thomson cooling profiles from HDF5 data.")
    
    parser.add_argument(
        "-p", "--property", 
        type=str, 
        default="pty", 
        help="Properties to plot. Include 'p' for pressure, 't' for temperature, 'y' for CO2 mass fraction. (Default: 'pty')"
    )
    
    parser.add_argument(
        "-r", "--reference", 
        action="store_true", 
        help="Flag to load and show reference data from 'odenburg.json'. (Default: False)"
    )

    args = parser.parse_args()

    #plot_property(args.property, args.reference)
