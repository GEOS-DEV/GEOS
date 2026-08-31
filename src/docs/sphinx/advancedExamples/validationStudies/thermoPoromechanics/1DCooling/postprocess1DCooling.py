"""LOCAL post-processing of the 1D confined cooling case (thermo-elastic vs Drucker-Prager).

Run by hand, after running both cases:

    geosx -i ThermoElastic_1DCooling_fim_smoke.xml        # in an elastic/ directory
    geosx -i ThermoDruckerPrager_1DCooling_fim_smoke.xml  # in a druckerPrager/ directory
    python3 postprocess1DCooling.py -e elastic -d druckerPrager -o <doc_directory>

It reads the GEOS TimeHistory outputs (stressHistory.hdf5, displacementHistory.hdf5) and
writes a single lightweight CSV:

    cooling1D.csv    temperature, GEOS stress and displacement of both models, analytical solution

That CSV is the file kept under version control. The .hdf5 files are not: the documentation
figure is produced by plot1DCooling.py, which reads this CSV only and depends solely on numpy
and matplotlib.

Analytical solution
-------------------
The column is fixed along y (both yneg and ypos are constrained) and free along x and z, hence
a state of uniaxial strain:
    eps_yy = 0,  sigma_xx = sigma_zz = 0
The thermo-elastic law then gives a stress that is purely induced by the cooling:
    sigma_yy = -E * alpha * dT          (tensile > 0 when dT < 0)

The Drucker-Prager criterion used by GEOS reads F = Q + b*P - c (DruckerPrager.hpp), with
P = tr(sigma)/3, Q the von Mises stress, and (DruckerPrager.cpp):
    b = 6 sin(phi) / (3 - sin(phi))
    c = 6 * cohesion * cos(phi) / (3 - sin(phi))
For this uniaxial state, Q = sigma_yy and P = sigma_yy/3, which yields the failure cap:
    sigma_failure = c / (1 + b/3)
and the cooling required to reach it: dT_failure = -sigma_failure / (E * alpha).
"""
import argparse
import math
import os
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from AnalyticalSol import DPAnalyticalSolution   # noqa: E402


# --- case properties (ThermoElastic/ThermoDruckerPrager_1DCooling_fim_smoke.xml) ---
BULK_MODULUS = 0.5e9        # defaultBulkModulus
SHEAR_MODULUS = 0.3e9       # defaultShearModulus
ALPHA = 3e-7                # defaultDrainedLinearTEC
COHESION = 5e3              # defaultCohesion
FRICTION_ANGLE = 15.27      # defaultFrictionAngle (degrees)
T_INITIAL = 100.0           # initialTemperature
T_FINAL = 20.0              # final value of timeFunction
T_RAMP_END = 100.0          # maxTime: end of the cooling ramp
COLUMN_WIDTH = 1.0          # extent along x (xCoords = { 0, 1 })


def youngModulus():
    K, G = BULK_MODULUS, SHEAR_MODULUS
    return 9.0 * K * G / (3.0 * K + G)


def druckerPragerCap():
    """(sigma_failure, dT_failure, b, c) for the uniaxial stress state of this case."""
    phi = math.radians(FRICTION_ANGLE)
    b = 6.0 * math.sin(phi) / (3.0 - math.sin(phi))
    c = 6.0 * COHESION * math.cos(phi) / (3.0 - math.sin(phi))
    sigma = c / (1.0 + b / 3.0)
    dT = -sigma / (youngModulus() * ALPHA)
    return sigma, dT, b, c


def imposedTemperature(t):
    """Ramp imposed by the TableFunction of the XML: 100 up to 1e-10, then linear down to 20."""
    return np.interp(t, [0.0, 1e-10, T_RAMP_END], [T_INITIAL, T_INITIAL, T_FINAL])


def elasticLateralDisplacement(dT):
    """Lateral displacement u_x of the free face, purely elastic branch.

    Same formulation as AnalyticalSol.DPAnalyticalSolution, without the plasticity check:
        eps_y^mech = -alpha*dT,   eps_x^mech = -eps_y^mech * lambda / (2 (lambda + G))
        eps_x      = eps_x^mech + alpha*dT
    The column is 1 m wide along x, hence u_x = eps_x * 1.
    """
    lam = BULK_MODULUS - 2.0 * SHEAR_MODULUS / 3.0
    epsMech_y = -ALPHA * dT
    epsMech_x = -epsMech_y * lam / (2.0 * (lam + SHEAR_MODULUS))
    return (epsMech_x + ALPHA * dT) * COLUMN_WIDTH


def _validLength(t):
    """The HDF5 buffer may be pre-allocated: cut at the first time value falling back to zero."""
    for j in range(1, len(t)):
        if t[j] < 1e-12:
            return j
    return len(t)


def readCase(case_dir):
    """(time, sigma_yy, u_x): stress at the centre of the column, lateral displacement of the
    free face x = 1. Both fields are uniform in this 1D setup, which is checked here."""
    ps = os.path.join(case_dir, "stressHistory.hdf5")
    pd = os.path.join(case_dir, "displacementHistory.hdf5")
    for p in (ps, pd):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p} not found. Run the corresponding case, then run this script again.")

    hs = h5py.File(ps, 'r')
    t = np.array(hs.get('averageStress Time')).ravel()
    S = np.array(hs.get('averageStress'))
    n = _validLength(t)
    t, S = t[:n], S[:n]
    spread = float(np.abs(S[-1, :, 1] - S[-1, :, 1].mean()).max())
    if spread > 1e-6 * max(1.0, abs(S[-1, :, 1].mean())):
        print(f"  [warning] {case_dir}: sigma_yy is not uniform (spread {spread:.3e} Pa)")
    sigma = S[:, S.shape[1] // 2, 1]

    hd = h5py.File(pd, 'r')
    td = np.array(hd.get('totalDisplacement Time')).ravel()[:n]
    U = np.array(hd.get('totalDisplacement'))[:n]
    X = np.array(hd.get('totalDisplacement ReferencePosition'))[0]
    if not np.allclose(t, td):
        raise ValueError(f"{case_dir}: stress and displacement time bases differ")
    free = np.where(np.abs(X[:, 0] - X[:, 0].max()) < 1e-9)[0]     # free face x = 1
    ux = U[:, free, 0].mean(axis=1)
    uy = float(np.abs(U[:, :, 1]).max())
    if uy > 1e-12:
        print(f"  [warning] {case_dir}: u_y is not zero ({uy:.3e} m), the column should be fixed")
    return t, sigma, ux


HEADER = ("time,temperature,dT,sigma_elastic_geos,sigma_druckerPrager_geos,"
          "sigma_elastic_analytical,sigma_yield,dp_b,dp_c,"
          "ux_elastic_geos,ux_druckerPrager_geos,"
          "ux_elastic_analytical,ux_druckerPrager_analytical")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-e", "--elastic", default="elastic", help="directory of the thermo-elastic run")
    ap.add_argument("-d", "--druckerPrager", default="druckerPrager", help="directory of the Drucker-Prager run")
    ap.add_argument("-o", "--outdir", default=".", help="directory where the CSV is written")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    te, se, uxe = readCase(args.elastic)
    td, sd, uxd = readCase(args.druckerPrager)
    if not np.allclose(te, td):
        raise ValueError("both runs do not share the same time base: "
                         "check that the Events of the two XML files are identical")

    T = imposedTemperature(te)
    dT = T - T_INITIAL
    E = youngModulus()
    sigma_ana = -E * ALPHA * dT                      # thermo-elastic, confined deformation
    sigma_yield, dT_yield, b, c = druckerPragerCap()

    # lateral displacement: elastic branch in closed form above, elasto-plastic branch from
    # AnalyticalSol.DPAnalyticalSolution (reference solution shipped with the PR)
    ux_ana_el = elasticLateralDisplacement(dT)
    dp = DPAnalyticalSolution()
    ux_ana_dp = np.array([dp.compute_disp(x)[0] * COLUMN_WIDTH for x in dT])

    # b and c are constant: stored as columns so that the plotting script can draw the failure
    # envelope without having to re-import the material properties
    data = np.column_stack([te, T, dT, se, sd, sigma_ana,
                            np.full_like(te, sigma_yield),
                            np.full_like(te, b), np.full_like(te, c),
                            uxe, uxd, ux_ana_el, ux_ana_dp])
    out = os.path.join(args.outdir, "cooling1D.csv")
    np.savetxt(out, data, delimiter=",", header=HEADER, comments="", fmt="%.8e")

    print(f"{out}  ({data.shape[0]} time steps)")
    print(f"  E = {E:.4g} Pa,  alpha = {ALPHA:.1e} 1/K")
    print(f"  Drucker-Prager: b = {b:.6f}, c = {c:.2f} Pa")
    print(f"  sigma_failure = {sigma_yield:.2f} Pa  reached at dT = {dT_yield:.2f} K "
          f"(T = {T_INITIAL + dT_yield:.2f})")
    err_el = np.abs(se - sigma_ana).max()
    print(f"  max deviation GEOS elastic / analytical: {err_el:.3e} Pa "
          f"({100 * err_el / max(abs(sigma_ana).max(), 1e-30):.3f} %)")
    print(f"  Drucker-Prager plateau, GEOS: {sd[-1]:.2f} Pa "
          f"(analytical {sigma_yield:.2f} Pa)")
    print(f"  final u_x, elastic       : GEOS {uxe[-1]:+.6e} m / analytical {ux_ana_el[-1]:+.6e} m "
          f"(deviation {abs(uxe[-1] - ux_ana_el[-1]):.2e})")
    print(f"  final u_x, Drucker-Prager: GEOS {uxd[-1]:+.6e} m / analytical {ux_ana_dp[-1]:+.6e} m "
          f"(deviation {abs(uxd[-1] - ux_ana_dp[-1]):.2e})")


if __name__ == "__main__":
    main()
