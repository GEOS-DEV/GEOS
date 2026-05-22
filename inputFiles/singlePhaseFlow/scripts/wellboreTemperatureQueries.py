"""
wellboreTemperatureQueries.py

Query script for the wellbore thermal diffusion validation tests.

Reads GEOS VTK output, extracts radial temperature profiles at selected timesteps, 
and writes model-results.txt that can be read by wellboreTemperatureFigure.py.

Usage:
    python wellboreTemperatureQueries.py --outputDir /path/to/vtkOutput [--snapshotIndices 1 2 5 10]

Output: model-results.txt with columns:
    time   r   temperature
"""

import os
import argparse

import numpy as np
from xml.etree import ElementTree

import vtk
from vtk.util.numpy_support import vtk_to_numpy


# ---------------------------------------------------------------------------
# VTK I/O helpers
# ---------------------------------------------------------------------------

def getTimestepsFromPVD(outputDir):
    """Parse vtkOutput.pvd → list of (time, vtmPath) tuples."""
    pvdPath = os.path.join(outputDir, 'vtkOutput.pvd')
    tree = ElementTree.parse(pvdPath)
    timesteps = []
    for dataset in tree.findall('.//DataSet'):
        time    = float(dataset.get('timestep'))
        vtmPath = os.path.join(outputDir, dataset.get('file'))
        timesteps.append((time, vtmPath))
    return timesteps


def _readVTU(vtuPath):
    """Read a single VTU file; return (centers, temperatures) as numpy arrays."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtuPath)
    reader.Update()
    grid = reader.GetOutput()
    cd   = grid.GetCellData()

    ec_vtk = cd.GetArray('elementCenter')
    if ec_vtk is None:
        raise RuntimeError(f"Array 'elementCenter' not found in {vtuPath}")
    centers = vtk_to_numpy(ec_vtk).reshape(-1, 3)

    t_vtk = cd.GetArray('temperature')
    if t_vtk is None:
        raise RuntimeError(f"Array 'temperature' not found in {vtuPath}")
    temperatures = vtk_to_numpy(t_vtk).ravel()

    return centers, temperatures


def readVTMData(vtmPath):
    """Read all VTU pieces in a VTM file; return one cell (closest to θ = 45°) per radial ring."""
    tree   = ElementTree.parse(vtmPath)
    vtmDir = os.path.dirname(vtmPath)

    all_centers, all_temps = [], []
    for ds in tree.findall('.//DataSet'):
        vtuFile = ds.get('file')
        if vtuFile is None:
            continue
        centers, temps = _readVTU(os.path.join(vtmDir, vtuFile))
        all_centers.append(centers)
        all_temps.append(temps)

    centers = np.vstack(all_centers)
    temps   = np.concatenate(all_temps)

    x, y         = centers[:, 0], centers[:, 1]
    r            = np.sqrt(x**2 + y**2)
    theta        = np.arctan2(y, x)
    target_theta = np.pi / 4.0

    r_rounded = np.round(r, decimals=8)
    unique_r  = np.unique(r_rounded)
    keep = []
    for rv in unique_r:
        mask = r_rounded == rv
        idx  = np.where(mask)[0]
        best = idx[np.argmin(np.abs(theta[idx] - target_theta))]
        keep.append(best)
    keep = np.array(keep)

    return r[keep], temps[keep]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Extract radial temperature profiles from GEOS VTK output.')
    parser.add_argument('--outputDir',
                        help='Path to GEOS VTK output directory (contains vtkOutput.pvd)',
                        default='.')
    parser.add_argument('--snapshotIndices', nargs='+', type=int,
                        default=[1, 2, 5, 10],
                        help='Indices into the PVD timestep list to extract (default: 1 2 5 10)')
    args = parser.parse_args()

    timesteps = getTimestepsFromPVD(args.outputDir)

    rows = []
    for idx in args.snapshotIndices:
        if idx >= len(timesteps):
            raise ValueError(f"Snapshot index {idx} out of range (file has {len(timesteps)} entries).")
        time, vtmPath = timesteps[idx]
        r, temperature = readVTMData(vtmPath)

        # Sort by radius
        sort_idx   = np.argsort(r)
        r          = r[sort_idx]
        temperature = temperature[sort_idx]

        for ri, Ti in zip(r, temperature):
            rows.append([time, ri, Ti])

    rows = np.array(rows)
    header = '         time             r   temperature'
    np.savetxt('model-results.txt', rows, header=header, fmt='%15.6g', comments='')
    print(f"Wrote model-results.txt ({len(rows)} rows, {len(args.snapshotIndices)} snapshots).")


if __name__ == '__main__':
    main()
