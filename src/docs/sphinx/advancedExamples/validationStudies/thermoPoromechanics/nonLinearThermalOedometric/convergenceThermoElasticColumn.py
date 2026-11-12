"""LOCAL mesh-convergence study of the column cooled by conduction (ThermoElasticColumnCooling_base.xml).

Run by hand; it runs GEOS itself on a sequence of meshes:

    python3 convergenceThermoElasticColumn.py -g <path to geosx> -w <work directory> -o <doc_directory>

For each number of cells nx it writes a deck that includes ThermoElasticColumnCooling_base.xml, with the time
step refined with the square of the cell size (backward Euler is first order in time, the finite volumes and
finite elements second order in space), runs it in <work directory>/nx<nx>, and compares the outputs with the
exact solution. It writes two lightweight CSV files, the ones kept under version control:

    columnCoolingConvergence.csv   cell size, time step and relative errors of each mesh
    columnCoolingProfile.csv       temperature and lateral stress along the column at the final time, coarse mesh

and exits with a non-zero status if an observed order falls below MIN_ORDER.

Exact solution
--------------
The fluid and the solid have the same volumetric heat capacity C and no fluid flows, so the energy balance is the
heat equation with diffusivity D = k / C, whatever the porosity. With T = T0 initially, T = T1 at x = 0 and an
insulated end at x = L:
    T(x, t) = T1 + ( T0 - T1 ) * sum_n 4 / ( (2n+1) pi ) * sin( l_n x ) * exp( -l_n^2 D t ),  l_n = (2n+1) pi / ( 2 L )
The column is confined in y and z and free at x = L, hence uniaxial strain with sigma_xx = 0 everywhere:
    sigma_yy = sigma_zz = - E / ( 1 - nu ) * I( x, t ),   eps_xx = 3 K / ( K + 4 G / 3 ) * I( x, t ),
with I( x, t ) the integral of alpha from T0 to T( x, t ). The end displacement is u_x( L ) = int_0^L eps_xx dx.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.abspath(os.path.join(HERE, "../../../../../../../inputFiles/thermoPoromechanics"))
BASE_DECK = os.path.join(INPUT_DIR, "ThermoElasticColumnCooling_base.xml")

# --- case properties (ThermoElasticColumnCooling_base.xml) ---
BULK_MODULUS = 9.12345679e9
SHEAR_MODULUS = 6.00813e9
CONDUCTIVITY = 1.0
HEAT_CAPACITY = 1.0e6            # volumetric, solid and fluid
T_INITIAL, T_COLD = 90.0, 0.0
LENGTH = 1.0
TABLE_KNOTS = np.array([0.0, 20.0, 50.0, 70.0])
TABLE_VALUES = np.array([0.2e-5, 0.5e-5, 0.6e-5, 1.6e-5])

# --- study ---
LEVELS = [10, 20, 40, 80, 160]   # number of cells along the column
DT_COARSE = 2500.0               # time step of the coarsest mesh; dt = DT_COARSE * ( LEVELS[0] / nx )^2
TIMES = [1.0e4, 2.5e4, 5.0e4]    # comparison times (multiples of every time step)
MIN_ORDER = 1.8

DECK = """<?xml version="1.0" ?>
<Problem>
  <Included>
    <File name="{base}"/>
  </Included>
  <Mesh>
    <InternalMesh name="mesh1" elementTypes="{{ C3D8 }}" xCoords="{{ 0, {length} }}" yCoords="{{ 0, 0.05 }}"
                  zCoords="{{ 0, 0.05 }}" nx="{{ {nx} }}" ny="{{ 1 }}" nz="{{ 1 }}" cellBlockNames="{{ column }}"/>
  </Mesh>
  <Events maxTime="{maxTime}">
    <PeriodicEvent name="temperatureHistoryCollection" timeFrequency="{every}" target="/Tasks/temperatureCollection"/>
    <PeriodicEvent name="stressHistoryCollection" timeFrequency="{every}" target="/Tasks/stressCollection"/>
    <PeriodicEvent name="displacementHistoryCollection" timeFrequency="{every}" target="/Tasks/displacementCollection"/>
    <PeriodicEvent name="historyOutput" timeFrequency="{every}" target="/Outputs/historyOutput"/>
    <PeriodicEvent name="displacementHistoryOutput" timeFrequency="{every}" target="/Outputs/displacementHistoryOutput"/>
    <!-- after the collections, so that each one records the state at its own time -->
    <PeriodicEvent name="solverApplications" forceDt="{dt}" target="/Solvers/thermoPoroSolver"/>
  </Events>
</Problem>
"""


def alpha(T):
    """Piecewise-linear table, constant outside its range."""
    return np.interp(T, TABLE_KNOTS, TABLE_VALUES)


def thermalStrain(T):
    """Integral of alpha from T_INITIAL to T, exact for the piecewise-linear table (trapezoids between knots)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    out = np.empty_like(T)
    for i, t in enumerate(T):
        lo, hi = sorted((T_INITIAL, t))
        x = np.concatenate(([lo], TABLE_KNOTS[(TABLE_KNOTS > lo) & (TABLE_KNOTS < hi)], [hi]))
        value = np.sum(0.5 * (alpha(x[1:]) + alpha(x[:-1])) * np.diff(x))
        out[i] = value if t >= T_INITIAL else -value
    return out


def exactTemperature(x, t, numTerms=400):
    diffusivity = CONDUCTIVITY / HEAT_CAPACITY
    n = np.arange(numTerms)[:, None]
    lam = (2 * n + 1) * np.pi / (2 * LENGTH)
    series = np.sum(4.0 / ((2 * n + 1) * np.pi) * np.sin(lam * x[None, :]) * np.exp(-lam ** 2 * diffusivity * t), axis=0)
    return T_COLD + (T_INITIAL - T_COLD) * series


def moduli():
    K, G = BULK_MODULUS, SHEAR_MODULUS
    E = 9 * K * G / (3 * K + G)
    nu = (3 * K - 2 * G) / (2 * (3 * K + G))
    return E, nu, 3 * K / (K + 4 * G / 3)


def exactEndDisplacement(t, numPoints=20001):
    _, _, strainFactor = moduli()
    x = np.linspace(0.0, LENGTH, numPoints)
    eps = strainFactor * thermalStrain(exactTemperature(x, t))
    return np.sum(0.5 * (eps[1:] + eps[:-1]) * np.diff(x))


def run(geos, workDir, nx):
    dt = DT_COARSE * (LEVELS[0] / nx) ** 2
    runDir = os.path.join(workDir, f"nx{nx}")
    os.makedirs(runDir, exist_ok=True)
    deck = os.path.join(runDir, "deck.xml")
    with open(deck, "w") as f:
        f.write(DECK.format(base=BASE_DECK, length=LENGTH, nx=nx, dt=dt, maxTime=TIMES[-1] + dt, every=DT_COARSE))
    with open(os.path.join(runDir, "log"), "w") as log:
        subprocess.run([geos, "-i", deck], cwd=runDir, stdout=log, stderr=subprocess.STDOUT, check=True)
    return runDir, dt


def read(runDir):
    with h5py.File(os.path.join(runDir, "columnHistory.hdf5"), "r") as h:
        time = np.array(h["temperature Time"])[:, 0]
        temperature = np.array(h["temperature"])
        stress = np.array(h["rockSolid_stress"])[:, :, 1]   # sigma_yy, first quadrature point
        centers = np.array(h["temperature elementCenter"])[0, :, 0]
    with h5py.File(os.path.join(runDir, "columnDisplacementHistory.hdf5"), "r") as h:
        displacement = np.array(h["totalDisplacement"])[:, :, 0]
        position = np.array(h["totalDisplacement ReferencePosition"])[0, :, 0]
    endDisplacement = displacement[:, np.isclose(position, LENGTH)].mean(axis=1)
    return time, centers, temperature, stress, endDisplacement


def errors(runDir):
    E, nu, _ = moduli()
    time, centers, temperature, stress, endDisplacement = read(runDir)
    result = {"temperature": 0.0, "stress": 0.0, "displacement": 0.0}
    for t in TIMES:
        i = int(np.argmin(np.abs(time - t)))
        assert abs(time[i] - t) < 1e-6 * t, f"no output at t = {t} in {runDir}"
        T = exactTemperature(centers, t)
        sigma = -E / (1 - nu) * thermalStrain(T)
        rms = lambda a: np.sqrt(np.mean(a ** 2))
        result["temperature"] = max(result["temperature"], rms(temperature[i] - T) / (T_INITIAL - T_COLD))
        result["stress"] = max(result["stress"], rms(stress[i] - sigma) / np.max(np.abs(sigma)))
        u = exactEndDisplacement(t)
        result["displacement"] = max(result["displacement"], abs(endDisplacement[i] - u) / abs(u))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--geos", required=True, help="path to the geosx executable")
    parser.add_argument("-w", "--work", required=True, help="directory where the runs are made")
    parser.add_argument("-o", "--output", default=HERE, help="directory where the CSV files are written")
    args = parser.parse_args()

    rows = []
    for nx in LEVELS:
        runDir, dt = run(os.path.abspath(args.geos), os.path.abspath(args.work), nx)
        e = errors(runDir)
        rows.append((LENGTH / nx, dt, e["temperature"], e["stress"], e["displacement"]))
        print(f"nx = {nx:4d}  dt = {dt:9.4f}  errors: temperature {e['temperature']:.3e}  stress {e['stress']:.3e}  "
              f"end displacement {e['displacement']:.3e}")
    rows = np.array(rows)
    np.savetxt(os.path.join(args.output, "columnCoolingConvergence.csv"), rows, delimiter=",", comments="",
               header="h,dt,error_temperature,error_stress,error_displacement", fmt="%.10e")

    ok = True
    for j, name in enumerate(["temperature", "stress", "end displacement"]):
        orders = np.log(rows[:-1, 2 + j] / rows[1:, 2 + j]) / np.log(rows[:-1, 0] / rows[1:, 0])
        fitted = np.polyfit(np.log(rows[:, 0]), np.log(rows[:, 2 + j]), 1)[0]
        print(f"observed order, {name:16s}: pairwise {np.array2string(orders, precision=2)}  fitted {fitted:.2f}")
        ok = ok and orders.min() >= MIN_ORDER

    # profile of the coarse run at the final time, for the figure
    time, centers, temperature, stress, _ = read(os.path.join(os.path.abspath(args.work), f"nx{LEVELS[1]}"))
    i = int(np.argmin(np.abs(time - TIMES[-1])))
    x = np.linspace(0.0, LENGTH, 201)
    E, nu, _ = moduli()
    T = exactTemperature(x, TIMES[-1])
    with open(os.path.join(args.output, "columnCoolingProfile.csv"), "w") as f:
        f.write("x,temperature_analytical,sigma_yy_analytical,x_geos,temperature_geos,sigma_yy_geos\n")
        for k in range(len(x)):
            geosPart = f"{centers[k]:.10e},{temperature[i, k]:.10e},{stress[i, k]:.10e}" if k < len(centers) else ",,"
            f.write(f"{x[k]:.10e},{T[k]:.10e},{-E / (1 - nu) * thermalStrain(T[k])[0]:.10e},{geosPart}\n")
    print("PASS" if ok else f"FAIL: an observed order is below {MIN_ORDER}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
