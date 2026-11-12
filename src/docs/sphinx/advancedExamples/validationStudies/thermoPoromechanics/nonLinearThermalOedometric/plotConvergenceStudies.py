"""Figures of the two mesh-convergence studies with a tabulated thermal expansion coefficient.

Executed by Sphinx (the `.. plot::` directives of Example.rst) on every documentation build.
It reads versioned CSV files sitting next to it:

    columnCoolingProfile.csv, columnCoolingConvergence.csv       (convergenceThermoElasticColumn.py)
    wellboreCoolingProfile.csv, wellboreCoolingConvergence.csv   (convergenceThermoElasticWellbore.py)

Those CSV files are produced offline by the two convergence scripts, which run GEOS. This script depends only on
numpy and matplotlib. The `.. plot::` directive calls the function column() or wellbore() named after the file;
run by hand, pass "column" or "wellbore" to draw one study, or nothing to draw both.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # the Sphinx `.. plot::` directive may execute this file without defining __file__
    HERE = os.getcwd()


def read(name):
    return np.genfromtxt(os.path.join(HERE, name), delimiter=",", names=True)


def slopeTwo(ax, h, error):
    """Reference line of slope 2 through the finest point."""
    ax.loglog(h, error[-1] * (h / h[-1]) ** 2, "k:", lw=1.5, label="Slope 2")


def column():
    profile, convergence = read("columnCoolingProfile.csv"), read("columnCoolingConvergence.csv")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    geos = ~np.isnan(profile["x_geos"])
    axes[0].plot(profile["x"], profile["temperature_analytical"], "b-", lw=2, label="Analytical")
    axes[0].plot(profile["x_geos"][geos], profile["temperature_geos"][geos], "k+", ms=10, mew=1.5, label="GEOS, 20 cells")
    axes[0].set_ylabel("Temperature")
    axes[1].plot(profile["x"], profile["sigma_yy_analytical"] * 1e-6, "b-", lw=2, label="Analytical")
    axes[1].plot(profile["x_geos"][geos], profile["sigma_yy_geos"][geos] * 1e-6, "k+", ms=10, mew=1.5, label="GEOS, 20 cells")
    axes[1].set_ylabel("Lateral stress $\\sigma_{yy}$ [MPa]")
    for ax in axes[:2]:
        ax.set_xlabel("x [m] (cold end at x = 0)")
        ax.set_title("Final time")
    h = convergence["h"]
    for key, marker, label in [("error_temperature", "o-", "Temperature"), ("error_stress", "s-", "Lateral stress"),
                               ("error_displacement", "^-", "End displacement")]:
        axes[2].loglog(h, convergence[key], "k" + marker, mfc="none", label=label)
    slopeTwo(axes[2], h, convergence["error_stress"])
    axes[2].set_xlabel("Cell size h [m] (time step proportional to $h^2$)")
    axes[2].set_ylabel("Relative error")
    axes[2].set_title("Mesh convergence")
    for ax in axes:
        ax.grid(True, which="both")
        ax.legend(loc="best")
    plt.tight_layout()


def wellbore():
    profile, convergence = read("wellboreCoolingProfile.csv"), read("wellboreCoolingConvergence.csv")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    geos = ~np.isnan(profile["r_geos"])
    axes[0].plot(profile["r"], profile["temperature_analytical"], "b-", lw=2, label="Analytical")
    axes[0].plot(profile["r_geos"][geos], profile["temperature_geos"][geos], "k+", ms=8, mew=1.2, label="GEOS, 8 x 8 cells")
    axes[0].set_ylabel("Temperature")
    for key, line, marker, label in [("sigma_tt", "b-", "k+", "$\\sigma_{\\theta\\theta}$"),
                                     ("sigma_zz", "b--", "kx", "$\\sigma_{zz}$"),
                                     ("sigma_rr", "b-.", "k.", "$\\sigma_{rr}$")]:
        axes[1].plot(profile["r"], profile[f"{key}_analytical"] * 1e-6, line, lw=2, label=f"{label}, analytical")
        axes[1].plot(profile["r_geos"][geos], profile[f"{key}_geos"][geos] * 1e-6, marker, ms=7, label=f"{label}, GEOS")
    axes[1].set_ylabel("Stress [MPa]")
    for ax in axes[:2]:
        ax.set_xlabel("Radius r [m] (wall at r = 0.1)")
        ax.set_title("Steady state")
    n = convergence["n"]
    h = 1.0 / n
    for key, marker, label in [("error_temperature", "o-", "Temperature"), ("error_sigma_rr", "s-", "$\\sigma_{rr}$"),
                               ("error_sigma_tt", "D-", "$\\sigma_{\\theta\\theta}$"),
                               ("error_wall_displacement", "^-", "Wall displacement")]:
        axes[2].loglog(h, convergence[key], "k" + marker, mfc="none", label=label)
    slopeTwo(axes[2], h, convergence["error_sigma_tt"])
    axes[2].set_xlabel("1 / number of cells in r and in $\\theta$")
    axes[2].set_ylabel("Relative error")
    axes[2].set_title("Mesh convergence")
    for ax in axes:
        ax.grid(True, which="both")
        ax.legend(loc="best", fontsize=9)
    plt.tight_layout()


def main():
    which = sys.argv[1:] or ["column", "wellbore"]
    if "column" in which:
        column()
    if "wellbore" in which:
        wellbore()
    plt.show()


if __name__ == "__main__":
    main()
