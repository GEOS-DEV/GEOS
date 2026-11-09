#!/usr/bin/env python3
"""Draw the load-displacement response of the single-edge notched benchmark under shear.

This is the light half of the post-processing: it reads two small text files sitting next
to it and needs only numpy and matplotlib, so the figure can be rebuilt with the
documentation, without the simulation output and without h5py.

  loadCurve.csv           GEOS result, produced from the TimeHistory output of the
                          benchmark deck by ../extractLoadCurve.py
  reference_Miehe2010.csv reference curve digitized from Miehe et al. (2010)
"""

import os

import numpy as np
import matplotlib.pyplot as plt

# When this script is executed by the "plot" directive of the documentation, Sphinx runs it
# without defining __file__, but with the working directory set to the directory holding
# it. Both situations are covered here.
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
GEOS_CURVE = os.path.join(HERE, "loadCurve.csv")
REFERENCE_CURVE = os.path.join(HERE, "reference_Miehe2010.csv")


def read_curve(path):
    """Read a two-column comma-separated curve, ignoring the commented provenance."""
    data = np.loadtxt(path, delimiter=",", comments="#")
    return data[:, 0], data[:, 1]


def main():
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)

    if os.path.isfile(REFERENCE_CURVE):
        displacement, load = read_curve(REFERENCE_CURVE)
        ax.plot(displacement, load, color="black", linewidth=1.8, label="Miehe et al. (2010)")

    displacement, load = read_curve(GEOS_CURVE)
    ax.plot(displacement,
            load,
            color="#cd1313",
            linestyle="none",
            marker="o",
            markersize=5,
            markeredgewidth=1.4,
            markerfacecolor="none",
            label="GEOS")

    ax.set_xlabel("Displacement $u$ [mm]", fontsize=12)
    ax.set_ylabel("Applied force [kN]", fontsize=12)
    ax.grid(True, color="#d9d9d9", linewidth=0.8)
    ax.legend(frameon=True, facecolor="white", edgecolor="black", fontsize=12, loc="lower right")
    plt.show()


if __name__ == "__main__":
    main()
