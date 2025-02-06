import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def read_hdf5(filename):
    """Reads an HDF5 file and prints its contents."""

    hf = h5py.File(filename, 'r')
    t = hf.get('phaseVolumeFraction Time')[:][:,0]
    x = hf.get('phaseVolumeFraction elementCenter')[:][0][:,0]
    s = hf.get('phaseVolumeFraction')[:][-1,:,0]
    hf.close()

    return t, x, s
 
if __name__ == "__main__":

    # Get the full path of the current working directory
    current_dir = os.getcwd()

    # Specify the folder name you want to access
    #folder_name = "my_folder"

    # Join the current directory path with the folder name
    #full_path = os.path.join(current_dir, folder_name)

    plt.figure(1)
    filename = os.path.join(current_dir, "results/saturationHistory.hdf5")  # Replace with the actual filename
    t1, x, s1 = read_hdf5(filename)
    plt.plot(x, s1, color='r', label = 'rho_c')
    print(t1)

    filename = os.path.join(current_dir, "resultsZ/saturationHistory.hdf5")  # Replace with the actual filename
    t2, x, s2 = read_hdf5(filename)
    plt.plot(x, s2, color='b', linestyle='--', label = 'z_c')
    print(t2)

    plt.xlabel("location")
    plt.ylabel("Gas saturation")
    plt.xlim(0, max(x))
    plt.ylim(0, 1)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15))
    plt.grid(True)
    plt.savefig('sat_profiles.png')

    plt.figure(2)
    plt.semilogy(x, s1-s2, color='g', label = 'z_c')
    plt.xlabel("location")
    plt.ylabel("saturation difference")
    plt.xlim(0, max(x))
    plt.grid(True)
    plt.savefig('sat_diff.png')