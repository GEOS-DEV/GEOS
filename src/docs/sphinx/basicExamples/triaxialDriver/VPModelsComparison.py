import os
import sys

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

def plotASingleModel(ax, path, label, color):
    # Load GEOSX results
    time, ax_strain, ra_strain1, ra_strain2, ax_stress, ra_stress1, ra_stress2, newton_iter, residual_norm = np.loadtxt(path, skiprows=5, unpack=True)

    p_num = -(ax_stress + 2.0 * ra_stress1) / 3.0 / 1.0e6
    q_num = -(ax_stress - ra_stress1) / 1.0e6
    strain_vol = ax_strain + 2.0 * ra_strain1

    #Visualization
    fsize = 30

    ax[0].plot(-ax_strain * 100,
               q_num,
               '-',
               color=color,
               mec=color,
               lineWidth=5)
    ax[0].plot(-ra_strain1 * 100, 
               q_num, 
               '-', 
               color=color,
               linewidth=5)

    ax[0].set_xlabel(r'Strain (%)', size=fsize, weight="bold")
    ax[0].set_ylabel(r'Deviatoric Stress (MPa)', size=fsize, weight="bold")
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)
    
    ax[1].plot(-ax_strain * 100,
               -strain_vol * 100,
               '-',
               color=color,
               linewidth=5,
               label=label)

    ax[1].set_xlabel(r'Axial Strain (%)', size=fsize, weight="bold")
    ax[1].set_ylabel(r'Volumetric Strain (%)', size=fsize, weight="bold")
    ax[1].legend(loc='lower center', bbox_to_anchor=(-0.3, 1.0),
          fancybox=True, shadow=True, ncol=5, fontsize=fsize)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

def main():
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(1, 2, figsize=(40, 10))

    plotASingleModel(ax, "DruckerPragerResults.txt", "DruckerPrager", cmap(0))
    plotASingleModel(ax, "ViscoDruckerPragerResults.txt",  "ViscoDruckerPrager", cmap(1))
    plotASingleModel(ax, "ExtendedDruckerPragerResults.txt", "ExtendedDruckerPrager", cmap(2))
    plotASingleModel(ax, "ViscoExtendedDruckerPragerResults.txt",  "ViscoExtendedDruckerPrager", cmap(3))

    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    plt.savefig("VPModelsComparion.png")
    os.system("xdg-open VPModelsComparion.png")

if __name__ == "__main__":
    main()
