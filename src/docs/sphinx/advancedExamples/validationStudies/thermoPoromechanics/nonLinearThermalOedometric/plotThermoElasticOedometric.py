"""Figure of the oedometric tests with a temperature-dependent thermal expansion coefficient.

Executed by Sphinx (the `.. plot::` directive of Example.rst) on every documentation build.
It reads a single versioned CSV file sitting next to it:

    thermoElasticOedometric.csv

That CSV is produced offline by postprocessThermoElasticOedometric.py from the TimeHistory (.hdf5)
outputs of the three benchmarks. This script therefore depends only on numpy and matplotlib.
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
    data = np.genfromtxt(os.path.join(HERE, "thermoElasticOedometric.csv"), delimiter=",", names=True)

    cases = [("linear", "Heating, linear $\\alpha(T)$"),
             ("tabulated", "Heating, tabulated $\\alpha(T)$"),
             ("cooling", "Cooling, linear $\\alpha(T)$")]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (name, title) in zip(axes, cases):
        temperature = data[f"temperature_{name}"]
        geos = data[f"sigma_xx_geos_{name}"]
        analytical = data[f"sigma_xx_analytical_{name}"]
        ax.plot(temperature, -analytical * 1e-5, "b-", lw=2, label="Analytical")
        ax.plot(temperature[::4], -geos[::4] * 1e-5, "k+", ms=8, mew=1.5, label="GEOS")
        ax.set_title(title)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Lateral compressive stress $-\\sigma_{xx}$ [bar]")
        ax.grid(True)
        ax.legend(loc="best")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
