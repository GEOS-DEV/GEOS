#!/usr/bin/env python3
"""
Figures of the "CO2 injection into a depleted reservoir" (dome) example.

Every ``plot_*`` function below creates one figure from the HDF5 file written
by ``dome_vtk_to_hdf5.py``.  They are called by the documentation through the
matplotlib ``plot`` directive of Sphinx, e.g.

    .. plot:: ../inputFiles/.../dome/plot_dome.py plot_pressure

so that no image is stored in the repository: the pictures (including the
rock-type map of the dome) are regenerated from the simulation output every
time the documentation is built.

Where things are read from
--------------------------
* results   : ``dome_results.hdf5`` next to this script, or the file given by
              the environment variable ``DOME_HDF5``;
* wells     : ``dome_mesh_bench.xml`` (``InternalWell`` polylines);
* schedule  : ``dome_soreide_whitson_base.xml`` (the three ``*.RATE`` tables)
              and the gas-water contact (``phaseContacts``).

Reading the wells, the schedule and the contact from the input files keeps the
figures consistent with the input if the example is modified.

Stand-alone use
---------------
    python plot_dome.py                       # write every figure as a PNG
    python plot_dome.py plot_pressure --hdf5 run/dome_results.hdf5 --outdir figs
"""
import argparse
import os
import re
import sys
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:                                    # __file__ is not always defined when
    _HERE = Path(__file__).resolve().parent   # run by the Sphinx plot directive
except NameError:                       # (the directive moves to the script dir)
    _HERE = Path.cwd()

SECONDS_PER_YEAR = 31557600.0

# --------------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------------
HDF5_NAME = "dome_results.hdf5"
MESH_XML = "dome_mesh_bench.xml"
PROPS_XML = "dome_properties.xml"
BASE_XML = "dome_soreide_whitson_base.xml"

# (time in years, caption).  The captions follow the schedule of the example:
# 0-8 yr production (+ water injection from 5 yr), 8-13 yr shut-in,
# 13-18 yr gas injection, 18-20 yr shut-in.
TIMES_ALL = [
    (0.0, "initial state"),
    (5.0, "end of depletion"),
    (8.0, "end of water flooding"),
    (13.0, "end of shut-in, start of gas injection"),
    (15.5, "2.5 yr into gas injection"),
    (18.0, "end of gas injection"),
    (20.0, "end of simulation"),
]
TIMES_CO2 = [
    (13.5, "0.5 yr into gas injection"),
    (15.5, "2.5 yr into gas injection"),
    (18.0, "end of gas injection"),
    (20.0, "end of simulation"),
]

ROCK_COLORS = ["#00aaff", "#ff0000", "#00ff00"]   # SAND, CLAY, CAPROCK
ROCK_NAMES = ["SAND", "CLAY", "CAPROCK"]
CAPROCK_GRAY = "#d9d9d9"

# Names of the datasets written by dome_vtk_to_hdf5.py (first match wins)
ALIASES = {
    "pressure": ["pressure"],
    "gas_saturation": ["phaseVolumeFraction_gas", "phaseVolumeFraction_0"],
    "CH4": ["globalCompFraction_CH4", "globalCompFraction_0"],
    "CO2": ["globalCompFraction_CO2", "globalCompFraction_1"],
}

# Fall-backs used only if the XML files cannot be found next to the script
DEFAULT_WELLS = {
    "WATER.INJECTOR": (125.0, 25.0, -4505.0),
    "GAS.INJECTOR": (125.0, 25.0, -4505.0),
    "PRODUCER": (2056.25, 25.0, -4126.25),
}
DEFAULT_GWC = -4400.0


# --------------------------------------------------------------------------
# Input-file helpers
# --------------------------------------------------------------------------
_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def _numbers(text):
    return np.array([float(v) for v in _NUM.findall(text)])


def read_wells(path=None):
    """Perforation mid-point (x, y, z) of every well of the mesh file."""
    path = Path(path) if path else _HERE / MESH_XML
    if not path.exists():
        return dict(DEFAULT_WELLS)
    wells = {}
    for w in ET.parse(path).getroot().iter("InternalWell"):
        pts = _numbers(w.attrib["polylineNodeCoords"]).reshape(-1, 3)
        perf = w.find("Perforation")
        d = float(perf.attrib["distanceFromHead"]) if perf is not None else 0.0
        v = pts[1] - pts[0]
        wells[w.attrib["wellRegionName"]] = tuple(pts[0] + v / np.linalg.norm(v) * d)
    return wells or dict(DEFAULT_WELLS)


def read_gwc(path=None):
    path = Path(path) if path else _HERE / BASE_XML
    if path.exists():
        for e in ET.parse(path).getroot().iter("HydrostaticEquilibrium"):
            return float(_numbers(e.attrib["phaseContacts"])[0])
    return DEFAULT_GWC


def read_table(name, path=None):
    """(coordinates, values) of a TableFunction of the base input file."""
    path = Path(path) if path else _HERE / BASE_XML
    for t in ET.parse(path).getroot().iter("TableFunction"):
        if t.attrib.get("name") == name:
            return _numbers(t.attrib["coordinates"]), _numbers(t.attrib["values"])
    raise KeyError(name)


def _rate_at(coords, vals, t):
    """Value of a table with ``interpolation="lower"`` (piecewise constant)."""
    i = np.searchsorted(coords, t, side="right") - 1
    return 0.0 if i < 0 else vals[i]


def read_schedule(path=None):
    """Return {'P','W','G': (coords_s, values)} for producer/water/gas."""
    try:
        return {"P": read_table("PRODUCER.RATE", path),
                "W": read_table("WATER.INJECTOR.RATE", path),
                "G": read_table("GAS.INJECTOR.RATE", path)}
    except (KeyError, OSError):
        return None


PHASE_LABELS = {frozenset(): "shut-in", frozenset("P"): "depletion",
                frozenset("PW"): "water flooding", frozenset("G"): "gas injection"}
PHASE_COLORS = {"shut-in": "#eeeeee", "depletion": "#fddbc7",
                "water flooding": "#d1e5f0", "gas injection": "#d9f0d3"}


def operating_phases(schedule, t_end):
    """[(t0, t1, label)] in seconds, deduced from the rate tables."""
    cuts = {0.0, t_end}
    for coords, _ in schedule.values():
        cuts.update(c for c in coords if 0.0 < c < t_end)
    cuts = sorted(cuts)
    out = []
    for t0, t1 in zip(cuts[:-1], cuts[1:]):
        tm = 0.5 * (t0 + t1)
        active = frozenset(k for k, (c, v) in schedule.items() if _rate_at(c, v, tm) > 0)
        out.append((t0, t1, PHASE_LABELS.get(active, "+".join(sorted(active)))))
    return out


# --------------------------------------------------------------------------
# Results access
# --------------------------------------------------------------------------
def hdf5_path():
    cands = [os.environ.get("DOME_HDF5"), _HERE / HDF5_NAME, Path.cwd() / HDF5_NAME]
    for c in cands:
        if c and Path(c).is_file():
            return Path(c)
    return None


class Results:
    """Thin reader that exposes the quasi-2D (x, z) section of the dome."""

    def __init__(self, path, y_target=None):
        import h5py
        self.h5 = h5py.File(path, "r")
        self.t_yr = self.h5["times_years"][:]
        m = self.h5["mesh"]
        self.center = m["center"][:]
        self.size = m["size"][:]
        self.region = m["region_id"][:]
        # section that contains the wells (the mesh is one cell thick in y in
        # the reference case; for several cells the one closest to y_target)
        ys = np.unique(self.center[:, 1].round(3))
        if y_target is None:
            y_target = ys[len(ys) // 2]
        self.y_sel = ys[np.argmin(np.abs(ys - y_target))]
        self.mask = np.abs(self.center[:, 1] - self.y_sel) < 1e-3
        cx = self.center[self.mask, 0].round(3)
        cz = self.center[self.mask, 2].round(3)
        self.xc, self.zc = np.unique(cx), np.unique(cz)
        self._ix = np.searchsorted(self.xc, cx)
        self._iz = np.searchsorted(self.zc, cz)
        self.xe, self.ze = _edges(self.xc), _edges(self.zc)

    def index(self, t_years):
        i = int(np.argmin(np.abs(self.t_yr - t_years)))
        if abs(self.t_yr[i] - t_years) > 0.3:
            warnings.warn(f"no output near {t_years} yr; using {self.t_yr[i]:.2f} yr")
        return i

    def dataset(self, key):
        for name in ALIASES.get(key, [key]):
            if name in self.h5["fields"]:
                return self.h5["fields"][name]
        raise KeyError(f"none of {ALIASES.get(key, [key])} found in the HDF5 file; "
                       "run dome_vtk_to_hdf5.py --list-fields on the VTK output")

    def to_grid(self, values_1d):
        g = np.full((len(self.zc), len(self.xc)), np.nan)
        g[self._iz, self._ix] = values_1d[self.mask]
        return g

    def field2d(self, key, it):
        return self.to_grid(self.dataset(key)[it, :].astype(float))

    def region2d(self):
        return self.to_grid(self.region.astype(float))


def _edges(c):
    if len(c) == 1:
        return np.array([c[0] - 0.5, c[0] + 0.5])
    e = np.empty(len(c) + 1)
    e[1:-1] = 0.5 * (c[1:] + c[:-1])
    e[0] = c[0] - 0.5 * (c[1] - c[0])
    e[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    return e


def _open():
    """Results object, or None (with a warning) if the file is missing."""
    p = hdf5_path()
    if p is None:
        warnings.warn(f"{HDF5_NAME} not found: run dome_vtk_to_hdf5.py on the GEOS "
                      "output or set DOME_HDF5")
        return None
    return Results(p, y_target=list(read_wells().values())[0][1])


def _placeholder(text="Results file not found"):
    fig, ax = plt.subplots(figsize=(8, 1.5))
    ax.text(0.5, 0.5, f"{text}\n(run dome_vtk_to_hdf5.py on the GEOS output "
            f"to create {HDF5_NAME})", ha="center", va="center", fontsize=9)
    ax.axis("off")
    return fig


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------
def _decorate(ax, R, wells, gwc, label=False):
    """Interfaces, gas-water contact and wells on one panel."""
    ax.set_facecolor(CAPROCK_GRAY)                    # inactive cells (CAPROCK)
    ax.contour(R.xc, R.zc, R.region2d(), levels=[0.5, 1.5], colors="k",
               linewidths=0.4, alpha=0.5)
    ax.axhline(gwc, color="k", ls="--", lw=0.7)
    inj = wells.get("GAS.INJECTOR", wells.get("WATER.INJECTOR"))
    prod = wells.get("PRODUCER")
    for w, marker, name in ((inj, "^", "injectors"), (prod, "v", "producer")):
        if w is not None:
            ax.plot(w[0], w[2], marker=marker, ms=6, mfc="k", mec="w", mew=0.8,
                    zorder=5, clip_on=False)
    if label:
        ax.annotate("gas-water contact", (R.xe[-1], gwc), xytext=(-4, 2),
                    textcoords="offset points", ha="right", va="bottom", fontsize=7)
        if inj is not None:
            ax.annotate("injectors", (inj[0], inj[2]), xytext=(6, 4),
                        textcoords="offset points", fontsize=7)
        if prod is not None:
            ax.annotate("producer", (prod[0], prod[2]), xytext=(-6, 4),
                        textcoords="offset points", fontsize=7, ha="right")
    ax.set_xlim(R.xe[0], R.xe[-1])
    ax.set_ylim(R.ze[0], R.ze[-1])
    ax.tick_params(labelsize=7)


def _fraction_limits(arrays):
    top = np.nanmax([np.nanmax(a) for a in arrays])
    top = 1.0 if not np.isfinite(top) or top <= 0 else min(1.0, np.ceil(top * 20) / 20)
    return 0.0, top


def _stack(key, times, cmap, cbar_label, diverging=False, as_change=False,
           fraction=False):
    """One figure: the key times stacked vertically, one shared colour scale."""
    R = _open()
    if R is None:
        return _placeholder()
    wells, gwc = read_wells(), read_gwc()
    its = [R.index(t) for t, _ in times]
    arrays = [R.field2d(key, i) for i in its]
    if key == "pressure":
        arrays = [a / 1e6 for a in arrays]                  # Pa -> MPa
    if as_change:
        ref = R.field2d(key, R.index(0.0)) / 1e6
        arrays = [a - ref for a in arrays]

    if diverging:
        m = np.nanmax([np.nanmax(np.abs(a)) for a in arrays]) or 1.0
        kw = dict(norm=TwoSlopeNorm(0.0, -m, m))
    elif fraction:
        vmin, vmax = _fraction_limits(arrays)
        kw = dict(vmin=vmin, vmax=vmax)
    else:
        kw = dict(vmin=np.nanmin([np.nanmin(a) for a in arrays]),
                  vmax=np.nanmax([np.nanmax(a) for a in arrays]))

    n = len(times)
    fig, axes = plt.subplots(n, 1, figsize=(10, 1.25 * n + 0.8), sharex=True,
                             layout="constrained")
    axes = np.atleast_1d(axes)
    for k, (ax, (t, cap), a, i) in enumerate(zip(axes, times, arrays, its)):
        pcm = ax.pcolormesh(R.xe, R.ze, a, cmap=cmap, shading="flat",
                            rasterized=True, **kw)
        _decorate(ax, R, wells, gwc, label=(k == 0))
        ax.set_title(f"t = {R.t_yr[i]:g} yr: {cap}", loc="left", fontsize=8, pad=2)
    axes[-1].set_xlabel("x (m)", fontsize=8)
    fig.supylabel("elevation z (m)", fontsize=8)
    fig.colorbar(pcm, ax=list(axes), label=cbar_label, pad=0.01, shrink=0.8, aspect=30)
    return fig


# --------------------------------------------------------------------------
# Public figures (one per ``.. plot::`` directive)
# --------------------------------------------------------------------------
def plot_rock_types():
    """Rock-type map of the dome, rebuilt from the cell regions of the output."""
    R = _open()
    if R is None:
        return _placeholder()
    wells, gwc = read_wells(), read_gwc()
    fig, ax = plt.subplots(figsize=(10, 2.9), layout="constrained")
    ax.pcolormesh(R.xe, R.ze, R.region2d(), cmap=ListedColormap(ROCK_COLORS),
                  vmin=-0.5, vmax=2.5, shading="flat", rasterized=True)
    _decorate(ax, R, wells, gwc, label=True)
    ax.set_xlabel("x (m)", fontsize=8)
    ax.set_ylabel("elevation z (m)", fontsize=8)
    ax.legend(handles=[Patch(color=c, label=n) for c, n in zip(ROCK_COLORS, ROCK_NAMES)]
              + [Line2D([], [], color="k", ls="--", lw=0.7, label="gas-water contact")],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4, frameon=False,
              fontsize=8)
    return fig


def plot_pressure():
    """Absolute pressure (MPa) at the key times."""
    return _stack("pressure", TIMES_ALL, "viridis", "pressure (MPa)")


def plot_pressure_change():
    """Pressure change with respect to the initial state (MPa)."""
    return _stack("pressure", TIMES_ALL[1:], "RdBu_r",
                  "pressure change from t = 0 (MPa)", diverging=True, as_change=True)


def plot_gas_saturation():
    """Gas saturation at the key times."""
    return _stack("gas_saturation", TIMES_ALL, "YlOrBr", "gas saturation (-)",
                  fraction=True)


def plot_co2_fraction():
    """CO2 global mole fraction during and after the gas injection."""
    return _stack("CO2", TIMES_CO2, "Greens", "CO$_2$ global mole fraction (-)",
                  fraction=True)


def plot_ch4_fraction():
    """CH4 global mole fraction at the key times."""
    return _stack("CH4", TIMES_ALL, "Purples", "CH$_4$ global mole fraction (-)",
                  fraction=True)


def plot_history():
    """Well-rate schedule (top) and pressure history (bottom)."""
    R = _open()
    if R is None:
        return _placeholder()
    sched = read_schedule()
    wells = read_wells()
    t_end = R.t_yr[-1] * SECONDS_PER_YEAR
    yrs = R.t_yr
    act = ~np.isnan(R.dataset("pressure")[0, :])
    p = R.dataset("pressure")[:, :].astype(float)[:, act] / 1e6
    vol = np.prod(R.size[act], axis=1)
    p_mean = (p * vol).sum(axis=1) / vol.sum()

    def cell_series(w):
        d = np.linalg.norm(R.center[act] - np.array(w), axis=1)
        return p[:, int(np.argmin(d))]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 5.4), sharex=True,
                                 height_ratios=[1, 1.5], layout="constrained")
    handles = []
    a1w = a1.twinx()                      # water rate is ~100x smaller: own axis
    if sched is not None:
        series = [("P", "producer (gas), left axis", "tab:red", a1),
                  ("G", "gas injector (85% CO$_2$, 15% CH$_4$), left axis", "tab:green", a1),
                  ("W", "water injector, right axis", "tab:blue", a1w)]
        for k, lab, col, ax in series:
            c, v = sched[k]
            xs = [0.0] + [x for x in c if 0.0 < x < t_end] + [t_end]
            ys = [_rate_at(c, v, x) for x in xs[:-1]]
            handles += ax.step(np.array(xs) / SECONDS_PER_YEAR, ys + [ys[-1]],
                               where="post", color=col, lw=1.8, label=lab)
        for (t0, t1, lab) in operating_phases(sched, t_end):
            for ax in (a1, a2):
                ax.axvspan(t0 / SECONDS_PER_YEAR, t1 / SECONDS_PER_YEAR,
                           color=PHASE_COLORS.get(lab, "#eeeeee"), alpha=0.6, lw=0)
            a1.text(0.5 * (t0 + t1) / SECONDS_PER_YEAR, 1.03, lab, ha="center",
                    va="bottom", fontsize=8, transform=a1.get_xaxis_transform())
        a1.legend(handles=handles, fontsize=7, loc="center right", frameon=True)
        a1.set_ylim(bottom=-0.05 * a1.get_ylim()[1])
        a1w.set_ylim(bottom=-0.05 * a1w.get_ylim()[1])
    a1w.set_ylabel("water rate (m$^3$/s)", fontsize=8, color="tab:blue")
    a1w.tick_params(labelsize=8, colors="tab:blue")
    a1.set_ylabel("target surface rate (m$^3$/s)", fontsize=8)
    a2.plot(yrs, p_mean, "k-", lw=2.2, label="volume-averaged, reservoir")
    if "GAS.INJECTOR" in wells:
        a2.plot(yrs, cell_series(wells["GAS.INJECTOR"]), color="tab:green", lw=1.2,
                label="injector cell")
    if "PRODUCER" in wells:
        a2.plot(yrs, cell_series(wells["PRODUCER"]), color="tab:red", lw=1.2,
                label="producer cell")
    a2.set_ylabel("pressure (MPa)", fontsize=8)
    a2.set_xlabel("time (years)", fontsize=8)
    a2.set_xlim(0, yrs[-1])
    a2.legend(fontsize=7, loc="best")
    for ax in (a1, a2):
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3, lw=0.4)
    return fig


def _scale(name, path=None):
    """``scale`` of a FieldSpecification of the properties file."""
    path = Path(path) if path else _HERE / PROPS_XML
    for e in ET.parse(path).getroot().iter("FieldSpecification"):
        if e.attrib.get("name") == name:
            return float(e.attrib["scale"])
    raise KeyError(name)


def plot_saturation_functions():
    """Relative permeabilities and capillary pressures of SAND and CLAY.

    Built from the tables of the input file, no simulation output needed.  The
    clay capillary curve is a J-function; it is converted to pressure with
    J * sigma * sqrt(porosity / permeability), using the surface tension of the
    model, the clay porosity and its permeability in the X direction.
    """
    sigma = None
    for e in ET.parse(_HERE / BASE_XML).getroot().iter("JFunctionCapillaryPressure"):
        sigma = float(e.attrib["wettingNonWettingSurfaceTension"])
    phi, kx = _scale("POROSITY.CLAY"), _scale("PERMX.CLAY")
    j_to_pa = sigma * np.sqrt(phi / kx)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    for rock, ls in (("SAND", "-"), ("CLAY", "--")):
        sw, krw = read_table(f"KRW.{rock}")
        sg, krg = read_table(f"KRG.{rock}")
        a1.plot(sw, krw, ls, color="tab:blue", label=f"water, {rock}")
        a1.plot(1.0 - sg, krg, ls, color="tab:red", label=f"gas, {rock}")
        sw_pc, pc = read_table(f"PCGW.{rock}")
        a2.plot(sw_pc, pc * (j_to_pa if rock == "CLAY" else 1.0) / 1e6, ls,
                color="k", label=rock)
    a1.set_xlabel("water saturation (-)")
    a1.set_ylabel("relative permeability (-)")
    a1.legend(fontsize=8)
    a2.set_xlabel("water saturation (-)")
    a2.set_ylabel("gas-water capillary pressure (MPa)")
    a2.legend(fontsize=8)
    for ax in (a1, a2):
        ax.set_xlim(0, 1)
        ax.grid(alpha=0.3, lw=0.4)
    return fig


ALL_FIGURES = ["plot_rock_types", "plot_saturation_functions", "plot_history", "plot_pressure",
               "plot_pressure_change", "plot_gas_saturation",
               "plot_co2_fraction", "plot_ch4_fraction"]


def _cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("figures", nargs="*", help=f"any of {ALL_FIGURES} (default: all)")
    ap.add_argument("--hdf5", help="results file (default: %s)" % HDF5_NAME)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--dpi", type=int, default=150)
    a = ap.parse_args(argv)
    if a.hdf5:
        os.environ["DOME_HDF5"] = a.hdf5
    os.makedirs(a.outdir, exist_ok=True)
    for name in a.figures or ALL_FIGURES:
        fig = globals()[name]()
        out = Path(a.outdir) / f"{name[5:]}.png"
        fig.savefig(out, dpi=a.dpi, bbox_inches="tight")
        plt.close(fig)
        print("wrote", out)


# The Sphinx plot directive runs this file with __name__ == "__main__", so the
# command-line entry point is disabled when Sphinx is loaded.
if __name__ == "__main__" and "sphinx" not in sys.modules:
    import matplotlib
    matplotlib.use("Agg")
    _cli()
