import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import sys
import argparse

# Import analytical functions from the top-level generate_analytical.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from generate_analytical import analytical_slip, analytical_traction, default_params


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # File path
    outputDir = args.outputDir
    hdf5File1Path = outputDir + "/displacementJump_ESG.hdf5"
    hdf5File2Path = outputDir + "/displacementJump_ISG.hdf5"

    # Read HDF5
    # Local Shear Displacement
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('displacementJump elementCenter')
    xl = np.array(xl)
    zcord_ESG = xl[0, :, 2]
    jump = hf.get('displacementJump')
    jump = np.array(jump)
    shearJump_ESG = np.abs(jump[-1, :, 1])
    hf.close()

    # Local Shear Displacement
    hf = h5py.File(hdf5File2Path, 'r')
    xl = hf.get('displacementJump elementCenter')
    xl = np.array(xl)
    zcord_ISG = xl[0, :, 2]
    jump = hf.get('displacementJump')
    jump = np.array(jump)
    shearJump_ISG = np.abs(jump[-1, :, 1])
    hf.close()

    # Compute analytical solution at a dense plotting grid
    params = default_params()
    z_ref = np.linspace(-300, 300, 2000)
    ts_ref = analytical_traction(z_ref, params) / 1.0e6   # Pa -> MPa
    ds_ref = analytical_slip(z_ref, params)

    fsize = 30
    msize = 12
    lw = 6
    malpha = 0.8
    lalpha = 1.0
    fig, ax = plt.subplots(1,figsize=(12, 16))
    cmap = plt.get_cmap("tab10")


    ax.plot(ds_ref, z_ref, linestyle='-', color=cmap(-1), alpha=lalpha, lw=lw, label='Analytical')
    ax.plot(shearJump_ESG, zcord_ESG, 'd--', color=cmap(1), markersize=msize, fillstyle='none', alpha=malpha, lw=lw, label='GEOS_ESG')
    ax.plot(shearJump_ISG, zcord_ISG, ':', color=cmap(2), markersize=msize, alpha=malpha, lw=lw, label='GEOS_ISG')
    ax.grid()
    ax.set_ylim(-300, 300)
    ax.set_xlabel('Shear Slip [m]', size=fsize, weight="bold")
    ax.set_ylabel('Depth [m]', size=fsize, weight="bold")
    ax.legend(loc='lower right', fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)

    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)

    plt.show()


if __name__ == "__main__":
    main()
