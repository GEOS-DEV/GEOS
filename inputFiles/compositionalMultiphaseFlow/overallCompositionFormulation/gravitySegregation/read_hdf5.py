import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np

def read_hdf5(filename):
    """Reads an HDF5 file and returns its contents."""
    with h5py.File(filename, 'r') as hf:
        t = hf.get('phaseVolumeFraction Time')[:][:,0]
        x = hf.get('phaseVolumeFraction elementCenter')[:][0][:,2]
        s = hf.get('phaseVolumeFraction')[:][:,:,0]  # Read all time snapshots
    return t, x, s

if __name__ == "__main__":
    current_dir = os.getcwd()

    filename = os.path.join(current_dir, "results/saturationHistory_drain.hdf5")
    t1, x, s1 = read_hdf5(filename)

    # Define the specific timesteps
    specific_times = [7200, 36000, 82800]

    # Convert specific times from seconds to hours and round off
    specific_times_hours = [round(time / 3600) for time in specific_times]

    # Find the indices of the specific timesteps
    selected_times = [time for time in specific_times if time in t1 and time in t2]

    # Define colors for different times
    red_shades = ['#FF0000', '#FF6666', '#FF9999']
    blue_shades = ['#0000FF', '#6666FF', '#9999FF']
    green_shades = ['#00FF00', '#66FF66', '#99FF99']

    plt.figure(1)
    for i, time in enumerate(selected_times):
        idx1 = np.where(t1 == time)[0][0]
        idx2 = np.where(t2 == time)[0][0]
        plt.plot(x, s2[idx2], color=blue_shades[i], linestyle='--', label=f'composition formulation t = {specific_times_hours[i]}h')

    plt.xlabel("Location")
    plt.ylabel("Gas Saturation")
    plt.xlim(0, max(x))
    plt.ylim(0, 1)
    plt.legend(loc='upper left', fontsize='x-small', bbox_to_anchor=(0.01, 0.85))
    plt.grid(True)
    plt.title("Gas Saturation Profiles")
    plt.savefig('sat_profiles.png')
