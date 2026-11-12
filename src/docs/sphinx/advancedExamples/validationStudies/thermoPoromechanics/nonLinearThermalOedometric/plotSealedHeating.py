"""Figure of the sealed oedometric heating with a tabulated thermal expansion coefficient.

Executed by Sphinx (the `.. plot::` directive of Example.rst) on every documentation build.
It reads a single versioned CSV file sitting next to it:

    sealedHeating.csv

That CSV is produced offline by postprocessSealedHeating.py from the TimeHistory (.hdf5) outputs of the two
smoke decks. This script therefore depends only on numpy and matplotlib.
"""
import os

import numpy as np
import matplotlib.pyplot as plt

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # the Sphinx `.. plot::` directive may execute this file without defining __file__
    HERE = os.getcwd()


def main():
    data = np.genfromtxt(os.path.join(HERE, "sealedHeating.csv"), delimiter=",", names=True)
    T = data["temperature"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    panels = [("pressure", 1e-6, "Pore pressure [MPa]"),
              ("porosity", 1.0, "Porosity [-]"),
              ("sigma_xx", 1e-6, "Lateral effective stress $\\sigma'_{xx}$ [MPa]")]
    for ax, (name, scale, label) in zip(axes, panels):
        ax.plot(T, data[f"{name}_analytical"] * scale, "b-", lw=2, label="Analytical")
        ax.plot(T, data[f"{name}_fim"] * scale, "k+", ms=10, mew=1.5, label="GEOS, fully implicit")
        ax.plot(T, data[f"{name}_sequential"] * scale, "ko", ms=8, mfc="none", label="GEOS, sequential")
        if name == "pressure":
            ax.plot(T, data["pressure_constantPorosityCoefficient"] * scale, "b--", lw=1.5,
                    label="Analytical, constant porosity coefficient")
        ax.set_xlabel("Temperature")
        ax.set_ylabel(label)
        ax.grid(True)
        ax.legend(loc="best")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
