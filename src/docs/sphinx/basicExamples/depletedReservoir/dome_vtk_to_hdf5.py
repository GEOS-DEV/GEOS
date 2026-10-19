#!/usr/bin/env python3
"""
Convert the GEOS VTK output of the dome example into a single HDF5 file.

GEOS writes one ``.vtm`` per output time (plus a ``.pvd`` index) and, inside
each ``.vtm``, one ``.vtu`` per MPI rank and per region.  This script reads that
tree and stores everything the documentation needs in one compact, rank
independent file:

    /times_s              (nt,)      output times [s]
    /times_years          (nt,)      output times [yr]  (1 yr = 31557600 s)
    /mesh/center          (nc, 3)    cell centres  x, y, z   [m]
    /mesh/size            (nc, 3)    cell extents dx, dy, dz [m]
    /mesh/region_id       (nc,)      0 = SAND, 1 = CLAY, 2 = CAPROCK ...
                                     (attribute ``region_names`` lists the names)
    /fields/<name>        (nt, nc)   one dataset per scalar field, NaN in the
                                     cells where the field does not exist
                                     (e.g. the inactive CAPROCK)

Field names
-----------
Scalar arrays keep their VTK name (``pressure``).  Multi-component arrays are
split into one dataset per component, named ``<array>_<component>``:

    phaseVolumeFraction_gas, phaseVolumeFraction_water
    globalCompFraction_CH4, globalCompFraction_CO2, globalCompFraction_H2S, ...

The component labels are taken from the VTK file when it carries them.  When
it does not, they are rebuilt from ``--phases`` / ``--components`` according
to the number of components of the array.

Usage
-----
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd -o dome_results.hdf5
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd --list-fields
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd --times-years 0 5 8 13 15.5 18 20
    python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd --fields pressure 'phaseVolumeFraction*'

Requirements: numpy, h5py, vtk (``pip install numpy h5py vtk``).
"""
import argparse
import datetime
import fnmatch
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

SECONDS_PER_YEAR = 31557600.0          # same Julian year as the XML time tables
DEFAULT_REGIONS = ("SAND", "CLAY", "CAPROCK")
DEFAULT_FIELDS = ("pressure", "temperature", "phaseVolumeFraction*",
                  "globalCompFraction*")
DEFAULT_PHASES = ("gas", "water")
DEFAULT_COMPONENTS = ("CH4", "CO2", "H2S", "H2O")


# --------------------------------------------------------------------------
# Index files
# --------------------------------------------------------------------------
def read_pvd(pvd_path):
    """Return a time-sorted list of (time_in_seconds, absolute_vtm_path)."""
    root = ET.parse(pvd_path).getroot()
    base = os.path.dirname(os.path.abspath(pvd_path))
    steps = []
    for ds in root.iter("DataSet"):
        steps.append((float(ds.attrib["timestep"]),
                      os.path.join(base, ds.attrib["file"])))
    steps.sort(key=lambda s: s[0])
    return steps


def region_files(vtm_path, regions):
    """Map region name -> list of .vtu paths (CellElementRegion blocks only)."""
    root = ET.parse(vtm_path).getroot()
    base = os.path.dirname(os.path.abspath(vtm_path))
    out = {}

    def walk(block, in_cell_regions, path):
        for child in block:
            if child.tag == "Block":
                name = child.attrib.get("name", "")
                walk(child, in_cell_regions or name == "CellElementRegion",
                     path + [name])
            elif child.tag == "DataSet" and in_cell_regions:
                # .../CellElementRegion/<REGION>/rank_i
                region = path[-1]
                if region in regions:
                    out.setdefault(region, []).append(
                        os.path.join(base, child.attrib["file"]))

    top = root.find("vtkMultiBlockDataSet")
    if top is None:
        raise RuntimeError(f"{vtm_path} is not a vtkMultiBlockDataSet file")
    walk(top, False, [])
    return out


# --------------------------------------------------------------------------
# VTK helpers (the vtk import is deferred so that --help works without it)
# --------------------------------------------------------------------------
def _vtk():
    try:
        import vtk
        from vtk.util import numpy_support
    except ImportError:
        sys.exit("This script needs the 'vtk' python package:  pip install vtk")
    return vtk, numpy_support


def read_vtu(path):
    vtk, _ = _vtk()
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def cell_centers(grid):
    vtk, ns = _vtk()
    f = vtk.vtkCellCenters()
    f.SetInputData(grid)
    f.Update()
    return ns.vtk_to_numpy(f.GetOutput().GetPoints().GetData()).astype(np.float64)


def cell_sizes(grid):
    """Axis-aligned extents of every cell (exact for the structured dome mesh)."""
    vtk, ns = _vtk()
    pts = ns.vtk_to_numpy(grid.GetPoints().GetData())
    n = grid.GetNumberOfCells()
    sizes = np.empty((n, 3))
    types = ns.vtk_to_numpy(grid.GetCellTypesArray())
    if np.all(types == vtk.VTK_HEXAHEDRON):
        conn = ns.vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(n, 8)
        p = pts[conn]                                   # (n, 8, 3)
        sizes[:] = p.max(axis=1) - p.min(axis=1)
    else:                                               # generic, slower
        for i in range(n):
            b = grid.GetCell(i).GetBounds()
            sizes[i] = (b[1] - b[0], b[3] - b[2], b[5] - b[4])
    return sizes


def component_labels(base, ncomp, vtk_array, phases, components):
    """Names of the components of a multi-component array."""
    labels = []
    for j in range(ncomp):
        name = vtk_array.GetComponentName(j)
        labels.append(name if name else None)
    if all(labels):
        return [l if l.startswith(base + "_") else f"{base}_{l}" for l in labels]
    if ncomp == len(phases):
        names = list(phases)
    elif ncomp == len(components):
        names = list(components)
    elif ncomp == len(phases) * len(components):       # [phase][component]
        names = [f"{p}_{c}" for p in phases for c in components]
    else:
        names = [str(j) for j in range(ncomp)]
    return [f"{base}_{n}" for n in names]


def read_region(files, phases, components, keep):
    """
    Read and concatenate all the rank files of one region.

    Returns (centers, sizes_or_None, {field_name: 1d array}) with the cells
    already sorted in a canonical (x, y, z) order, independent of the
    domain decomposition.
    """
    vtk, ns = _vtk()
    chunks = []
    for fn in sorted(files):
        g = read_vtu(fn)
        n = g.GetNumberOfCells()
        if n == 0:
            continue
        sel = np.ones(n, dtype=bool)                    # drop ghost cells
        cd = g.GetCellData()
        for ghost_name, owned in (("vtkGhostType", lambda a: a == 0),
                                  ("ghostRank", lambda a: a < 0)):
            arr = cd.GetArray(ghost_name)
            if arr is not None:
                sel &= owned(ns.vtk_to_numpy(arr).reshape(n, -1)[:, 0])
        rec = {"centers": cell_centers(g)[sel], "sizes": cell_sizes(g)[sel],
               "fields": {}}
        for i in range(cd.GetNumberOfArrays()):
            arr = cd.GetAbstractArray(i)
            if not isinstance(arr, vtk.vtkDataArray):
                continue
            base = arr.GetName()
            if not any(fnmatch.fnmatch(base, pat) for pat in keep):
                continue
            data = ns.vtk_to_numpy(arr).astype(np.float64).reshape(n, -1)[sel]
            if data.shape[1] == 1:
                rec["fields"][base] = data[:, 0]
            else:
                labs = component_labels(base, data.shape[1], arr, phases, components)
                for j, lab in enumerate(labs):
                    rec["fields"][lab] = data[:, j]
        chunks.append(rec)

    if not chunks:
        return np.empty((0, 3)), np.empty((0, 3)), {}

    centers = np.concatenate([c["centers"] for c in chunks])
    sizes = np.concatenate([c["sizes"] for c in chunks])
    names = set().union(*[c["fields"].keys() for c in chunks])
    fields = {k: np.concatenate([c["fields"].get(k, np.full(len(c["centers"]), np.nan))
                                 for c in chunks]) for k in names}

    # canonical order: sort by rounded (x, y, z)
    r = np.round(centers, 3)
    order = np.lexsort((r[:, 2], r[:, 1], r[:, 0]))
    return centers[order], sizes[order], {k: v[order] for k, v in fields.items()}


def read_step(vtm_path, regions, phases, components, keep):
    """Read one output time; returns the concatenation over regions."""
    files = region_files(vtm_path, regions)
    missing = [r for r in regions if r not in files]
    if missing:
        raise RuntimeError(f"{vtm_path}: region(s) {missing} not found")
    cen, siz, rid, flds = [], [], [], []
    for k, name in enumerate(regions):
        c, s, f = read_region(files[name], phases, components, keep)
        cen.append(c)
        siz.append(s)
        rid.append(np.full(len(c), k, dtype=np.int8))
        flds.append(f)
    counts = [len(c) for c in cen]
    names = sorted(set().union(*[f.keys() for f in flds]))
    fields = {}
    for nm in names:
        fields[nm] = np.concatenate(
            [f.get(nm, np.full(n, np.nan)) for f, n in zip(flds, counts)])
    return (np.concatenate(cen), np.concatenate(siz), np.concatenate(rid), fields)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def select_times(steps, years):
    if not years:
        return steps
    t = np.array([s[0] for s in steps]) / SECONDS_PER_YEAR
    picked = []
    for y in years:
        i = int(np.argmin(np.abs(t - y)))
        if abs(t[i] - y) > 0.26:                       # outputs are every 0.5 yr
            print(f"warning: no output within 0.25 yr of {y} yr "
                  f"(closest: {t[i]:.2f} yr)", file=sys.stderr)
        if steps[i] not in picked:
            picked.append(steps[i])
    return sorted(picked, key=lambda s: s[0])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pvd", help="GEOS VTK index file (e.g. VTK.OUTPUT.pvd)")
    ap.add_argument("-o", "--output", default="dome_results.hdf5")
    ap.add_argument("--regions", nargs="+", default=list(DEFAULT_REGIONS),
                    help="cell regions to export (default: %(default)s)")
    ap.add_argument("--fields", nargs="+", default=list(DEFAULT_FIELDS),
                    help="array-name patterns to export; use '*' for everything "
                         "(default: %(default)s)")
    ap.add_argument("--phases", nargs="+", default=list(DEFAULT_PHASES))
    ap.add_argument("--components", nargs="+", default=list(DEFAULT_COMPONENTS))
    ap.add_argument("--times-years", nargs="+", type=float, default=None,
                    help="only export the outputs closest to these times [yr]")
    ap.add_argument("--list-fields", action="store_true",
                    help="print the arrays found in the first output and exit")
    args = ap.parse_args(argv)

    steps = select_times(read_pvd(args.pvd), args.times_years)
    if not steps:
        sys.exit("no time step found in " + args.pvd)

    keep = ["*"] if args.list_fields else args.fields
    print(f"reading {len(steps)} time steps from {args.pvd}")
    cen0, siz0, rid0, f0 = read_step(steps[0][1], args.regions, args.phases,
                                     args.components, keep)
    if args.list_fields:
        print(f"{len(cen0)} cells; fields in the first output:")
        for k in sorted(f0):
            n_ok = int(np.sum(~np.isnan(f0[k])))
            print(f"  {k:45s} defined in {n_ok} cells")
        return 0

    import h5py
    names = sorted(f0)
    if not names:
        sys.exit("none of the requested fields was found; try --list-fields")
    nt, nc = len(steps), len(cen0)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["title"] = "GEOS dome example - CO2 injection into a depleted reservoir"
        h5.attrs["source_pvd"] = os.path.abspath(args.pvd)
        h5.attrs["created"] = datetime.datetime.now().isoformat(timespec="seconds")
        h5.attrs["phases"] = list(args.phases)
        h5.attrs["components"] = list(args.components)
        h5.attrs["seconds_per_year"] = SECONDS_PER_YEAR
        times = np.array([s[0] for s in steps])
        h5.create_dataset("times_s", data=times)
        h5.create_dataset("times_years", data=times / SECONDS_PER_YEAR)
        g = h5.create_group("mesh")
        g.create_dataset("center", data=cen0, compression="gzip")
        g.create_dataset("size", data=siz0, compression="gzip")
        g.create_dataset("region_id", data=rid0, compression="gzip")
        g.attrs["region_names"] = list(args.regions)
        fg = h5.create_group("fields")
        dsets = {k: fg.create_dataset(k, shape=(nt, nc), dtype="f4",
                                      chunks=(1, nc), compression="gzip",
                                      compression_opts=4, shuffle=True,
                                      fillvalue=np.nan)
                 for k in names}
        for k in names:
            dsets[k][0, :] = f0[k]

        for it in range(1, nt):
            print(f"  [{it + 1:3d}/{nt}] t = {steps[it][0] / SECONDS_PER_YEAR:7.2f} yr",
                  flush=True)
            cen, _, rid, flds = read_step(steps[it][1], args.regions, args.phases,
                                          args.components, keep)
            if len(cen) != nc or not np.allclose(cen, cen0, atol=1e-6):
                sys.exit(f"cell layout of {steps[it][1]} differs from the first "
                         "output; the mesh must be fixed in time")
            for k in names:
                if k in flds:
                    dsets[k][it, :] = flds[k]
    print(f"wrote {args.output}: {nt} times, {nc} cells, {len(names)} fields")
    return 0


if __name__ == "__main__":
    sys.exit(main())
