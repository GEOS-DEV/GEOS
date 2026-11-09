#!/usr/bin/env python3
"""Reduce the GEOS TimeHistory output of the single-edge notched shear benchmark to a
light load curve.

This is the heavy half of the post-processing. Run it, without any argument, from this
directory once the benchmark has been run and its history file copied here:

    python3 extractLoadCurveShear.py

It reads the history file directly where GEOS wrote it, under
``inputFiles/phaseField/benchmark/singleEdgeNotch``, and writes ``loadCurve.csv``, a
two-column text file of a few kilobytes holding the imposed
displacement and the resulting load. The figure is then drawn from that file alone by
``plotLoadCurve.py``, which needs neither h5py nor the simulation output.

The reaction force is the traction integrated over the loaded edge. The area of each
boundary face is its cell volume divided by the height of the cell row, itself recovered
from the distance between the edge and the cell centers, so no dimension of the mesh is
duplicated here. The load is reported per millimetre of thickness, as in the reference.
"""

import argparse
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# The history file is read where GEOS writes it, in the directory of the benchmark deck:
# seven levels above this script is the root of the GEOS repository.
GEOS_DIR = os.path.normpath(os.path.join(HERE, *[os.pardir] * 7))
DECK_DIR = os.path.join(GEOS_DIR, "inputFiles", "phaseField", "benchmark", "singleEdgeNotch")

HISTORY_FILE = "singleEdgeNotchShear_topBoundary.hdf5"
OUTPUT_FILE = "loadCurve.csv"

# Candidate locations, in order: the output directory of the run (geosx -o Shear),
# then the directory of the deck itself when the run was launched without -o.
HISTORY_CANDIDATES = [os.path.join(DECK_DIR, "Shear", HISTORY_FILE),
                      os.path.join(DECK_DIR, HISTORY_FILE),
                      os.path.join(HERE, HISTORY_FILE)]

# Voigt ordering of averageStress is (XX, YY, ZZ, YZ, XZ, XY). On a face of normal y, the
# shear test is driven by index 5 (sigma_xy), and the imposed displacement is its
# component 0.
STRESS_COMPONENT = 5
DISPLACEMENT_COMPONENT = 0

# The reference reports the load carried by a 1 mm thick specimen.
REFERENCE_THICKNESS_MM = 1.0


def find_history():
    """Return the first candidate history file that exists."""
    for candidate in HISTORY_CANDIDATES:
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError("{} not found in any of:\n  ".format(HISTORY_FILE)
                            + "\n  ".join(HISTORY_CANDIDATES))


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
    """Return the imposed displacement [mm] and the load [kN per mm of thickness]."""
    with h5py.File(path, "r") as history:
        time = dataset(history, "averageStress", "Time")
        stress = dataset(history, "averageStress")
        volume = dataset(history, "elementVolume")
        center = dataset(history, "averageStress", "elementCenter")
        position = dataset(history, "totalDisplacement", "ReferencePosition")[0]
        displacement = dataset(history, "totalDisplacement")

    # GEOS preallocates the datasets for the whole run, so an interrupted run leaves
    # unwritten samples behind. Keep only the rows whose collected time is increasing.
    stalled = np.flatnonzero(np.diff(time[:, 0]) <= 0.0)
    written = int(stalled[0]) + 1 if stalled.size else len(time)
    stress, displacement = stress[:written], displacement[:written]

    row_height = 2.0 * (position[:, 1].max() - center[0, :, 1].mean())
    face_area = volume[0] / row_height
    thickness = position[:, 2].max() - position[:, 2].min()

    force = (stress[:, :, STRESS_COMPONENT] * face_area).sum(axis=1)
    imposed = displacement[:, :, DISPLACEMENT_COMPONENT].mean(axis=1)

    load = force / thickness * REFERENCE_THICKNESS_MM * 1.0e-3
    return imposed, load


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
    imposed, load = read_history(history)
    imposed, load = imposed[::args.subsample], load[::args.subsample]

    header = ("Load curve of the single-edge notched benchmark under shear\n"
              "Extracted from {} by extractLoadCurveShear.py\n"
              "displacement [mm],load [kN per mm of thickness]").format(
                  os.path.basename(history))
    np.savetxt(args.out, np.column_stack((imposed, load)), delimiter=",",
               header=header, comments="# ", fmt="%.6e")

    peak = int(np.argmax(load))
    print("{} samples written to {}".format(len(load), args.out))
    print("peak load {:.4f} kN at a displacement of {:.4e} mm".format(load[peak], imposed[peak]))


if __name__ == "__main__":
    main()
