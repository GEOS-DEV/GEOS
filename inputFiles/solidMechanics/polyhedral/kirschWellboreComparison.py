"""
Kirsch wellbore: finite elements and mixed virtual elements against the analytical solution.

Companion of kirschWellbore_fem.xml and kirschWellbore_mixedVEM.xml. The analytical solution is
imported from the validation script of the documentation, so both comparisons use the same formulas.

These decks load the far field as a traction on the outer boundary instead of an initial stress,
so stresses are total and compare directly, while the displacement also carries the uniform
plane strain field of the far field. That field is subtracted, leaving the displacement induced by
the hole and the well pressure, which is what the analytical solution describes.

Usage, from the directory holding the output folders:
    python3 kirschWellboreComparison.py --fem kirsch_fem --mixedVEM kirsch_vem
"""

import argparse
import glob
import importlib.util
import os
import sys
import types
import xml.etree.ElementTree as ElementTree

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

scriptDir = os.path.dirname(os.path.abspath(__file__))
geosDir = os.path.abspath(os.path.join(scriptDir, "..", "..", ".."))
figureScript = os.path.join(geosDir, "src", "docs", "sphinx", "advancedExamples", "validationStudies",
                            "wellboreProblems", "kirschWellbore", "kirschWellboreFigure.py")

METHODS = ("FEM", "mixed VEM")
MARKERS = {"FEM": "o", "mixed VEM": "s"}
STRESS_NAMES = [r"$\sigma_{rr}$", r"$\sigma_{\theta\theta}$", r"$\sigma_{r\theta}$"]
DISP_NAMES = [r"$u_r$", r"$u_\theta$"]


def loadDocumentationScript():
    # the analytical class does not read HDF5, so a missing h5py must not block the import
    # a stub from an earlier call has no __spec__, which find_spec rejects
    if "h5py" not in sys.modules and importlib.util.find_spec("h5py") is None:
        sys.modules["h5py"] = types.ModuleType("h5py")
    spec = importlib.util.spec_from_file_location("kirschWellboreFigure", figureScript)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def readParameters(baseXml, loadXml):
    base = ElementTree.parse(baseXml)
    material = base.find("Constitutive/ElasticIsotropic")
    mesh = base.find("Mesh/InternalWellbore")
    radius = [float(v) for v in mesh.get("radius").strip("{} ").split(",")]

    params = {"mechanical": {"bulkModulus": float(material.get("defaultBulkModulus")),
                             "shearModulus": float(material.get("defaultShearModulus"))},
              "rw": radius[0],
              "rout": radius[1],
              "plotFileRoot": base.find("Outputs/VTK").get("plotFileRoot", "vtkOutput"),
              "stress": None,
              "wellLoad": None}

    for traction in ElementTree.parse(loadXml).findall("FieldSpecifications/Traction"):
        if traction.get("tractionType") == "stress":
            params["stress"] = np.array([float(v) for v in traction.get("inputStress").strip("{} ").split(",")])
        elif traction.get("tractionType") == "normal":
            params["wellLoad"] = float(traction.get("scale"))
    if params["stress"] is None or params["wellLoad"] is None:
        sys.exit(f"{loadXml}: expected a stress traction and a normal traction")

    return params


def readLastStep(outputDir, plotFileRoot):
    root = os.path.join(outputDir, plotFileRoot)
    files = sorted(glob.glob(os.path.join(root, "*", "**", "*.vtu"), recursive=True))
    if not files:
        sys.exit(f"no VTU output under {root}")
    step = lambda f: os.path.relpath(f, root).split(os.sep)[0]
    lastStep = max(step(f) for f in files)
    cells, points = {}, {}

    for f in files:
        if step(f) != lastStep:
            continue
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(f)
        reader.Update()
        grid = reader.GetOutput()

        # ghosts repeat the cells and nodes owned by another rank
        for store, data, coords in ((cells, grid.GetCellData(), None),
                                    (points, grid.GetPointData(), vtk_to_numpy(grid.GetPoints().GetData()))):
            owned = vtk_to_numpy(data.GetArray("ghostRank")) < 0
            if coords is not None:
                store.setdefault("coordinates", []).append(coords[owned])
            for i in range(data.GetNumberOfArrays()):
                store.setdefault(data.GetArrayName(i), []).append(vtk_to_numpy(data.GetArray(i))[owned])

    return ({k: np.concatenate(v) for k, v in cells.items()},
            {k: np.concatenate(v) for k, v in points.items()})


def polarStress(voigt, theta, rotate):
    # voigt order (xx, yy, zz, yz, xz, xy), the order both solvers write
    s = voigt
    tensor = np.array([[s[0], s[5], s[4]], [s[5], s[1], s[3]], [s[4], s[3], s[2]]])
    local = rotate(tensor, theta)
    return local[0][0], local[1][1], local[0][1]


def farFieldDisplacement(x, stress, mechanical):
    K, G = mechanical["bulkModulus"], mechanical["shearModulus"]
    nu = (3 * K - 2 * G) / (2 * (3 * K + G))
    # plane strain, tension positive
    exx = ((1 - nu) * stress[0] - nu * stress[1]) / (2 * G)
    eyy = ((1 - nu) * stress[1] - nu * stress[0]) / (2 * G)
    return np.array([exx * x[0], eyy * x[1], 0.0])


def binned(r, values, edges):
    # tetrahedral cell centers scatter around the profile, so average them in radial bins
    index = np.digitize(r, edges) - 1
    keep = [b for b in range(len(edges) - 1) if np.any(index == b)]
    return (np.array([np.mean(r[index == b]) for b in keep]),
            np.array([values[index == b].mean(axis=0) for b in keep]))


def exactSolution(documentation, params, theta0):
    # the documentation class takes compressive stresses and the well pressure as positive numbers
    analytical = documentation.Analytical(params["mechanical"], params["rw"], -params["stress"],
                                          -params["wellLoad"], theta0)

    def exact(r, theta):
        analytical.theta = theta
        # back to tension positive, the sign convention of the solvers
        return (-analytical.computeRadialStress(r), -analytical.computeHoopStress(r),
                -analytical.computeShearStress(r), analytical.computeRadialDisp(r), analytical.computeShearDisp(r))

    return exact


def extractProfile(documentation, params, exact, outputDir, theta0, band):
    cells, points = readLastStep(outputDir, params["plotFileRoot"])
    stressKey = "averageStress" if "averageStress" in cells else "stress"

    centers = cells["elementCenter"]
    cellTheta = np.arctan2(centers[:, 1], centers[:, 0])
    selected = np.abs(cellTheta - theta0) < np.radians(band)

    rS, numS, exS = [], [], []
    for i in np.flatnonzero(selected):
        r = np.hypot(centers[i, 0], centers[i, 1])
        rS.append(r)
        numS.append(polarStress(cells[stressKey][i], cellTheta[i], documentation.stressRotation))
        exS.append(exact(r, cellTheta[i])[:3])

    # nodal displacement for finite elements, the cell value u_h(x_E) for the mixed form
    if "totalDisplacement" in points:
        x, u = points["coordinates"], points["totalDisplacement"]
        t = np.arctan2(x[:, 1], x[:, 0])
        top = np.isclose(x[:, 2], x[:, 2].max())
        # the node ray closest to theta0, since theta0 need not be a multiple of the angular spacing
        nearest = t[top][np.argmin(np.abs(t[top] - theta0))]
        onProfile = top & (np.abs(t - nearest) < 1e-6)
    else:
        x, u, t = centers, cells["displacement"], cellTheta
        onProfile = selected

    rU, numU, exU = [], [], []
    for i in np.flatnonzero(onProfile):
        r = np.hypot(x[i, 0], x[i, 1])
        induced = u[i] - farFieldDisplacement(x[i], params["stress"], params["mechanical"])
        local = documentation.dispRotation(induced, t[i])
        rU.append(r)
        numU.append((local[0], local[1]))
        exU.append(exact(r, t[i])[3:])

    if not rS or not rU:
        sys.exit(f"{outputDir}: empty profile at theta = {np.degrees(theta0):g} deg, widen --band")

    profile = dict(zip(("rS", "numS", "exS", "rU", "numU", "exU"),
                       map(np.array, (rS, numS, exS, rU, numU, exU))))
    profile["errors"] = [np.linalg.norm(profile["numS"][:, c] - profile["exS"][:, c]) /
                         np.linalg.norm(profile["exS"][:, c]) for c in range(3)] + \
                        [np.linalg.norm(profile["numU"][:, c] - profile["exU"][:, c]) /
                         np.linalg.norm(profile["exU"][:, c]) for c in range(2)]
    return profile


def compareRuns(runs, params, theta, band, bins, save, show=False, title=None):
    """Plot one row per method, stress left and displacement right, and return the errors per method."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    documentation = loadDocumentationScript()
    theta0 = np.radians(theta)
    exact = exactSolution(documentation, params, theta0)

    rLine = np.geomspace(params["rw"], params["rout"], 400)
    curve = np.array([exact(r, theta0) for r in rLine])
    edges = np.geomspace(params["rw"], params["rout"], bins + 1)

    fsize, msize = 18, 9
    cmap = plt.get_cmap("tab10")
    fig, axes = plt.subplots(len(runs), 2, figsize=(22, 8 * len(runs)), sharex=True, sharey="col", squeeze=False)

    print(f"profile at theta = {theta} deg, cells within +/- {band} deg")
    print(f"{'method':12s} {'points':>7s} {'s_rr':>9s} {'s_tt':>9s} {'s_rt':>9s} {'u_r':>9s} {'u_t':>9s}   relative L2 error")

    allErrors = {}
    for row, (label, outputDir) in enumerate(runs):
        profile = extractProfile(documentation, params, exact, outputDir, theta0, band)
        errors = profile["errors"]
        allErrors[label] = errors
        print(f"{label:12s} {len(profile['rS']):7d} " + " ".join(f"{e:9.2e}" for e in errors))

        rSb, numSb = binned(profile["rS"], profile["numS"], edges)
        rUb, numUb = binned(profile["rU"], profile["numU"], edges)
        axStress, axDisp = axes[row]

        for c in range(3):
            axStress.semilogx(rLine, curve[:, c] / 1e6, lw=3, alpha=0.6, color=cmap(c),
                              label=STRESS_NAMES[c] + " analytical")
            axStress.semilogx(rSb, numSb[:, c] / 1e6, MARKERS[label], ms=msize, mfc="none", mew=1.8, color=cmap(c),
                              label=STRESS_NAMES[c] + " " + label)
        for c in range(2):
            axDisp.semilogx(rLine, curve[:, 3 + c] * 1e3, lw=3, alpha=0.6, color=cmap(c),
                            label=DISP_NAMES[c] + " analytical")
            axDisp.semilogx(rUb, numUb[:, c] * 1e3, MARKERS[label], ms=msize, mfc="none", mew=1.8, color=cmap(c),
                            label=DISP_NAMES[c] + " " + label)

        axStress.set_title(label + "   relative L2 error: " +
                           ", ".join(f"{n} {e:.1e}" for n, e in zip(STRESS_NAMES, errors[:3])), size=fsize * 0.85)
        axDisp.set_title(label + "   relative L2 error: " +
                         ", ".join(f"{n} {e:.1e}" for n, e in zip(DISP_NAMES, errors[3:])), size=fsize * 0.85)
        axStress.set_ylabel(r"$\sigma$ (MPa)", size=fsize)
        axDisp.set_ylabel("induced displacement (mm)", size=fsize)

    for a in axes.ravel():
        a.set_xlim(params["rw"], params["rout"])
        a.grid(True, which="both", alpha=0.3)
        a.tick_params(labelsize=fsize * 0.8)
        a.legend(fontsize=fsize * 0.6, ncol=2, loc="lower right")
    for a in axes[-1]:
        a.set_xlabel("r (m)", size=fsize)

    fig.suptitle(title or rf"Kirsch wellbore, $\theta$ = {theta:g}$^\circ$", size=fsize * 1.1)
    fig.tight_layout()
    fig.savefig(save, dpi=110)
    print(f"figure written to {save}")
    if show:
        plt.show()
    plt.close(fig)
    return allErrors


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fem", help="output directory of kirschWellbore_fem.xml")
    parser.add_argument("--mixedVEM", help="output directory of kirschWellbore_mixedVEM.xml")
    parser.add_argument("--base", default=os.path.join(scriptDir, "kirschWellbore_base.xml"),
                        help="deck holding the material, the mesh and the output")
    parser.add_argument("--loads", default=os.path.join(scriptDir, "kirschWellbore_fem.xml"),
                        help="deck holding the far field and well tractions")
    parser.add_argument("--theta", type=float, default=45.0, help="angle of the profile, degrees")
    parser.add_argument("--band", type=float, default=0.45, help="half width of the cell selection, degrees")
    parser.add_argument("--bins", type=int, default=30, help="log spaced radial bins for the markers")
    parser.add_argument("--save", default="kirschWellboreComparison.png", help="figure file")
    parser.add_argument("--show", action="store_true", help="also open the figure")
    args = parser.parse_args()

    runs = [(label, path) for label, path in zip(METHODS, (args.fem, args.mixedVEM)) if path]
    if not runs:
        parser.error("give at least one of --fem, --mixedVEM")

    compareRuns(runs, readParameters(args.base, args.loads), args.theta, args.band, args.bins, args.save, args.show)


if __name__ == "__main__":
    main()
