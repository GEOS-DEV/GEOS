"""LOCAL mesh-convergence study of the steady cooling of a wellbore wall (ThermoElasticWellboreSteadyCooling_base.xml).

Run by hand; it runs GEOS itself on a sequence of meshes:

    python3 convergenceThermoElasticWellbore.py -g <path to geosx> -w <work directory> -o <doc_directory>

For each level it writes a deck that includes ThermoElasticWellboreSteadyCooling_base.xml with nr = nt = n cells,
runs it in <work directory>/n<n>, and compares the steady state with the exact solution. It writes two lightweight
CSV files, the ones kept under version control:

    wellboreCoolingConvergence.csv   number of cells and relative errors of each mesh
    wellboreCoolingProfile.csv       temperature and stresses against the radius, exact and coarse mesh

and exits with a non-zero status if an observed order falls below MIN_ORDER.

Exact solution
--------------
Steady conduction in the annulus a < r < b with T( a ) = T1 and T( b ) = T0 (reached after a few long steps):
    T( r ) = T1 + ( T0 - T1 ) ln( r / a ) / ln( b / a )
The thermal strain e( r ) is the integral of alpha from T0 to T( r ). In plane strain, with lambda and G the Lame
coefficients and beta = 3 K / ( lambda + 2 G ), the radial displacement and the stresses are
    u( r )       = beta J( r ) / r + C1 r + C2 / r,          J( r ) = int_a^r e( s ) s ds
    sigma_rr     = -2 G beta J / r^2 + 2 ( lambda + G ) C1 - 2 G C2 / r^2
    sigma_tt     = -2 G beta e + 2 G beta J / r^2 + 2 ( lambda + G ) C1 + 2 G C2 / r^2
    sigma_zz     = lambda ( beta e + 2 C1 ) - 3 K e
with C1 and C2 given by the free inner wall, sigma_rr( a ) = 0, and the fixed outer boundary, u( b ) = 0.
J is computed by quadrature on a fine grid that contains the radii where T crosses the knots of the table.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.abspath(os.path.join(HERE, "../../../../../../../inputFiles/thermoPoromechanics"))
BASE_DECK = os.path.join(INPUT_DIR, "ThermoElasticWellboreSteadyCooling_base.xml")

# --- case properties (ThermoElasticWellboreSteadyCooling_base.xml) ---
BULK_MODULUS = 9.12345679e9
SHEAR_MODULUS = 6.00813e9
T_FAR, T_WALL = 90.0, 0.0
INNER_RADIUS, OUTER_RADIUS = 0.1, 1.0
TABLE_KNOTS = np.array([0.0, 20.0, 50.0, 70.0])
TABLE_VALUES = np.array([0.2e-5, 0.5e-5, 0.6e-5, 1.6e-5])

# --- study ---
LEVELS = [8, 16, 32, 64]   # nr = nt
MIN_ORDER = 1.8

DECK = """<?xml version="1.0" ?>
<Problem>
  <Included>
    <File name="{base}"/>
  </Included>
  <Mesh>
    <InternalWellbore name="mesh1" elementTypes="{{ C3D8 }}" radius="{{ {a}, {b} }}" theta="{{ 0, 90 }}"
                      zCoords="{{ 0, 0.01 }}" nr="{{ {n} }}" nt="{{ {n} }}" nz="{{ 1 }}"
                      trajectory="{{ {{ 0.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.01 }} }}" autoSpaceRadialElems="{{ 1 }}"
                      cellBlockNames="{{ rock }}"/>
  </Mesh>
  <Events maxTime="6e9">
    <PeriodicEvent name="solverApplications" forceDt="1e9" target="/Solvers/thermoPoroSolver"/>
    <PeriodicEvent name="historyCollection" timeFrequency="1e9" target="/Tasks/temperatureCollection"/>
    <PeriodicEvent name="stressHistoryCollection" timeFrequency="1e9" target="/Tasks/stressCollection"/>
    <PeriodicEvent name="displacementHistoryCollection" timeFrequency="1e9" target="/Tasks/displacementCollection"/>
    <PeriodicEvent name="historyOutput" timeFrequency="1e9" target="/Outputs/historyOutput"/>
    <PeriodicEvent name="displacementHistoryOutput" timeFrequency="1e9" target="/Outputs/displacementHistoryOutput"/>
  </Events>
</Problem>
"""


def alpha(T):
    return np.interp(T, TABLE_KNOTS, TABLE_VALUES)


def thermalStrain(T):
    """Integral of alpha from T_FAR to T, exact for the piecewise-linear table."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    out = np.empty_like(T)
    for i, t in enumerate(T):
        lo, hi = sorted((T_FAR, t))
        x = np.concatenate(([lo], TABLE_KNOTS[(TABLE_KNOTS > lo) & (TABLE_KNOTS < hi)], [hi]))
        value = np.sum(0.5 * (alpha(x[1:]) + alpha(x[:-1])) * np.diff(x))
        out[i] = value if t >= T_FAR else -value
    return out


def exactTemperature(r):
    return T_WALL + (T_FAR - T_WALL) * np.log(r / INNER_RADIUS) / np.log(OUTER_RADIUS / INNER_RADIUS)


class ExactSolution:
    """Plane-strain thermo-elastic annulus with a radial thermal strain profile."""

    def __init__(self, numPoints=400001):
        K, G = BULK_MODULUS, SHEAR_MODULUS
        self.K, self.G, self.lam = K, G, K - 2 * G / 3
        self.beta = 3 * K / (self.lam + 2 * G)
        a, b = INNER_RADIUS, OUTER_RADIUS
        kinkRadii = a * np.exp(np.log(b / a) * (TABLE_KNOTS - T_WALL) / (T_FAR - T_WALL))
        r = np.union1d(np.linspace(a, b, numPoints), kinkRadii[(kinkRadii > a) & (kinkRadii < b)])
        f = thermalStrain(exactTemperature(r)) * r
        self.r = r
        self.J = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))))
        Jb = self.J[-1]
        matrix = np.array([[2 * (self.lam + G), -2 * G / a ** 2], [b, 1.0 / b]])
        self.C1, self.C2 = np.linalg.solve(matrix, [0.0, -self.beta * Jb / b])

    def fields(self, r):
        G, lam, beta = self.G, self.lam, self.beta
        e = thermalStrain(exactTemperature(r))
        J = np.interp(r, self.r, self.J)
        u = beta * J / r + self.C1 * r + self.C2 / r
        srr = -2 * G * beta * J / r ** 2 + 2 * (lam + G) * self.C1 - 2 * G * self.C2 / r ** 2
        stt = -2 * G * beta * e + 2 * G * beta * J / r ** 2 + 2 * (lam + G) * self.C1 + 2 * G * self.C2 / r ** 2
        szz = lam * (beta * e + 2 * self.C1) - 3 * self.K * e
        return u, srr, stt, szz


def run(geos, workDir, n):
    runDir = os.path.join(workDir, f"n{n}")
    os.makedirs(runDir, exist_ok=True)
    deck = os.path.join(runDir, "deck.xml")
    with open(deck, "w") as f:
        f.write(DECK.format(base=BASE_DECK, a=INNER_RADIUS, b=OUTER_RADIUS, n=n))
    with open(os.path.join(runDir, "log"), "w") as log:
        subprocess.run([geos, "-i", deck], cwd=runDir, stdout=log, stderr=subprocess.STDOUT, check=True)
    return runDir


def read(runDir):
    """Final steady state: cell radius, temperature, polar stresses (mean of the quadrature points), wall displacement."""
    with h5py.File(os.path.join(runDir, "wellboreHistory.hdf5"), "r") as h:
        temperature = np.array(h["temperature"])[-1]
        centers = np.array(h["temperature elementCenter"])[-1]
        stress = np.array(h["rockSolid_stress"])[-1].reshape(len(temperature), -1, 6).mean(axis=1)
    with h5py.File(os.path.join(runDir, "wellboreDisplacementHistory.hdf5"), "r") as h:
        displacement = np.array(h["totalDisplacement"])[-1]
        position = np.array(h["totalDisplacement ReferencePosition"])[-1]
    r = np.hypot(centers[:, 0], centers[:, 1])
    c, s = centers[:, 0] / r, centers[:, 1] / r
    sxx, syy, szz, sxy = stress[:, 0], stress[:, 1], stress[:, 2], stress[:, 5]
    srr = c * c * sxx + s * s * syy + 2 * c * s * sxy
    stt = s * s * sxx + c * c * syy - 2 * c * s * sxy
    nodeRadius = np.hypot(position[:, 0], position[:, 1])
    wall = np.isclose(nodeRadius, INNER_RADIUS)
    wallDisplacement = np.mean((displacement[wall, 0] * position[wall, 0] + displacement[wall, 1] * position[wall, 1])
                               / nodeRadius[wall])
    return r, temperature, srr, stt, szz, wallDisplacement


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--geos", required=True, help="path to the geosx executable")
    parser.add_argument("-w", "--work", required=True, help="directory where the runs are made")
    parser.add_argument("-o", "--output", default=HERE, help="directory where the CSV files are written")
    args = parser.parse_args()
    exact = ExactSolution()
    rms = lambda x: np.sqrt(np.mean(x ** 2))

    rows = []
    for n in LEVELS:
        r, T, srr, stt, szz, uWall = read(run(os.path.abspath(args.geos), os.path.abspath(args.work), n))
        u, srrExact, sttExact, szzExact = exact.fields(r)
        uWallExact = exact.fields(np.array([INNER_RADIUS]))[0][0]
        scale = np.max(np.abs(sttExact))
        row = (n, rms(T - exactTemperature(r)) / (T_FAR - T_WALL), rms(srr - srrExact) / scale,
               rms(stt - sttExact) / scale, rms(szz - szzExact) / scale, abs(uWall - uWallExact) / abs(uWallExact))
        rows.append(row)
        print(f"n = {n:3d}  errors: temperature {row[1]:.3e}  sigma_rr {row[2]:.3e}  sigma_tt {row[3]:.3e}  "
              f"sigma_zz {row[4]:.3e}  wall displacement {row[5]:.3e}")
    rows = np.array(rows)
    np.savetxt(os.path.join(args.output, "wellboreCoolingConvergence.csv"), rows, delimiter=",", comments="",
               header="n,error_temperature,error_sigma_rr,error_sigma_tt,error_sigma_zz,error_wall_displacement",
               fmt="%.10e")

    ok = True
    for j, name in enumerate(["temperature", "sigma_rr", "sigma_tt", "sigma_zz", "wall displacement"]):
        orders = np.log(rows[:-1, 1 + j] / rows[1:, 1 + j]) / np.log(rows[1:, 0] / rows[:-1, 0])
        fitted = -np.polyfit(np.log(rows[:, 0]), np.log(rows[:, 1 + j]), 1)[0]
        print(f"observed order, {name:17s}: pairwise {np.array2string(orders, precision=2)}  fitted {fitted:.2f}")
        ok = ok and orders.min() >= MIN_ORDER

    # profile of the coarse run, for the figure
    r, T, srr, stt, szz, _ = read(os.path.join(os.path.abspath(args.work), f"n{LEVELS[0]}"))
    rFine = np.geomspace(INNER_RADIUS, OUTER_RADIUS, 201)
    _, srrExact, sttExact, szzExact = exact.fields(rFine)
    order = np.argsort(r)
    with open(os.path.join(args.output, "wellboreCoolingProfile.csv"), "w") as f:
        f.write("r,temperature_analytical,sigma_rr_analytical,sigma_tt_analytical,sigma_zz_analytical,"
                "r_geos,temperature_geos,sigma_rr_geos,sigma_tt_geos,sigma_zz_geos\n")
        for k in range(max(len(rFine), len(r))):
            exactPart = (f"{rFine[k]:.10e},{exactTemperature(rFine[k]):.10e},{srrExact[k]:.10e},{sttExact[k]:.10e},"
                         f"{szzExact[k]:.10e}") if k < len(rFine) else ",,,,"
            i = order[k] if k < len(r) else None
            geosPart = f"{r[i]:.10e},{T[i]:.10e},{srr[i]:.10e},{stt[i]:.10e},{szz[i]:.10e}" if i is not None else ",,,,"
            f.write(exactPart + "," + geosPart + "\n")
    print("PASS" if ok else f"FAIL: an observed order is below {MIN_ORDER}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
