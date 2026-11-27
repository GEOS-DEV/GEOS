#!/usr/bin/env python3
"""Reduce the GEOS TimeHistory output of the three-point bending with holes benchmark to a light load curve.

This is the heavy half of the post-processing. Run it, without any argument, once the
benchmark has been run:

    python3 extractLoadCurveThreePointsBendingWithHoles.py

It reads the history file directly where GEOS wrote it, under
``inputFiles/phaseField/benchmark/threePointsBendingWithHoles``, and writes ``loadCurve.csv``, a two-column
text file of a few kilobytes holding the deflection and the applied load. The figure is
then drawn from that file alone by ``plotLoadCurve.py``, which needs neither h5py nor the
simulation output.

The load is applied on a line of nodes at mid-span, so it cannot be recovered by
integrating the traction over a boundary face. It is obtained instead from the shear force
on a vertical cut between the left support and the load line, where the shear force of a
simply supported beam under a central load equals half of the applied load:

    P = 2 * A_cut * < sigma_xy >_V

with < . >_V the volume-weighted average over the cells of the cut. The area of the cut is
recovered from the history file itself: the load line lies on the top surface of a beam
whose lower face is at y = 0, so its position gives the height of the section, and the
extent of its reference positions gives the thickness.
"""

import argparse
import glob
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Seven levels above this script is the root of the GEOS repository.
GEOS_DIR = os.path.normpath(os.path.join(HERE, *[os.pardir] * 7))
DECK_DIR = os.path.join(GEOS_DIR, "inputFiles", "phaseField", "benchmark", "threePointsBendingWithHoles")

HISTORY_FILE = "threePointsBendWithHoles_reactionCut.hdf5"
OUTPUT_FILE = "loadCurve.csv"

# Voigt ordering of averageStress is (XX, YY, ZZ, YZ, XZ, XY): index 5 is sigma_xy, the
# vertical traction component on a cut of normal x.
SHEAR_COMPONENT = 5
DEFLECTION_COMPONENT = 1

# The reference reports the load carried by a 1 mm thick specimen.
REFERENCE_THICKNESS_MM = 1.0


def find_history():
    """Return the first existing candidate: the directory of the deck, an output directory
    created there by ``geosx -o``, or this directory."""
    candidates = ([os.path.join(DECK_DIR, HISTORY_FILE)]
                  + sorted(glob.glob(os.path.join(DECK_DIR, "*", HISTORY_FILE)))
                  + [os.path.join(HERE, HISTORY_FILE)])
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError("{} not found in any of:\n  ".format(HISTORY_FILE)
                            + "\n  ".join(candidates))


def dataset(history, field, kind=None):
    """Datasets are named "<field> [<meta>] <setName>", so they are looked up by field and
    kind rather than by set name, which makes this independent of the deck's set names."""
    if kind == "Time":
        return history[field + " Time"][:]
    for key in history:
        if not key.startswith(field + " ") or key.endswith(" Time"):
            continue
        meta = "elementCenter" in key or "ReferencePosition" in key
        if (kind is None) != meta:
            return history[key][:]
    raise KeyError("no {} dataset for field {}".format(kind or "value", field))


def read_history(path):
    """Return the deflection [mm] and the applied load [kN]."""
    with h5py.File(path, "r") as history:
        time = dataset(history, "averageStress", "Time")
        stress = dataset(history, "averageStress")
        volume = dataset(history, "elementVolume")
        position = dataset(history, "totalDisplacement", "ReferencePosition")[0]
        displacement = dataset(history, "totalDisplacement")

    # GEOS preallocates the datasets for the whole run, so a run that is still going or was
    # interrupted leaves unwritten samples behind. Keep the rows whose time is increasing.
    stalled = np.flatnonzero(np.diff(time[:, 0]) <= 0.0)
    written = int(stalled[0]) + 1 if stalled.size else len(time)
    stress, volume, displacement = stress[:written], volume[:written], displacement[:written]

    height = position[:, 1].max()
    thickness = position[:, 2].max() - position[:, 2].min()
    cut_area = height * thickness

    mean_shear = (stress[:, :, SHEAR_COMPONENT] * volume).sum(axis=1) / volume.sum(axis=1)
    load = 2.0 * cut_area * np.abs(mean_shear)
    deflection = np.abs(displacement[:, :, DEFLECTION_COMPONENT].mean(axis=1))

    # Reported per millimetre of thickness and in kN, as in the reference and as in the
    # single-edge notched benchmarks, so that the four cases are directly comparable.
    load = load / thickness * REFERENCE_THICKNESS_MM * 1.0e-3
    return deflection, load


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--history", default=None,
                        help="GEOS TimeHistory file (default: the run output of the deck).")
    parser.add_argument("--out", default=os.path.join(HERE, OUTPUT_FILE),
                        help="Load curve to write (default: %(default)s).")
    parser.add_argument("--subsample", type=int, default=1, help="Keep every n-th sample.")
    args = parser.parse_args()

    history = args.history if args.history else find_history()
    deflection, load = read_history(history)
    deflection, load = deflection[::args.subsample], load[::args.subsample]

    header = ("Load-deflection curve of the three-point bending with holes benchmark\n"
              "Extracted from {} by extractLoadCurveThreePointsBendingWithHoles.py\n"
              "deflection [mm],load [kN per mm of thickness]").format(os.path.basename(history))
    np.savetxt(args.out, np.column_stack((deflection, load)), delimiter=",",
               header=header, comments="# ", fmt="%.6e")

    peak = int(np.argmax(load))
    print("{} samples written to {}".format(len(load), args.out))
    print("peak load {:.4f} kN at a deflection of {:.4e} mm".format(load[peak], deflection[peak]))


if __name__ == "__main__":
    main()
