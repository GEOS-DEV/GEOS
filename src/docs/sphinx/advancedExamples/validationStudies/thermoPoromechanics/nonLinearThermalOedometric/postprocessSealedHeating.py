"""LOCAL post-processing of the sealed oedometric heating (ThermoPoroElasticSealedHeating_*_smoke.xml).

Run by hand, after running both decks, each in its own directory:

    geos -i ThermoPoroElasticSealedHeating_fim_smoke.xml          # in a fim/ directory
    geos -i ThermoPoroElasticSealedHeating_sequential_smoke.xml   # in a sequential/ directory
    python3 postprocessSealedHeating.py -f fim -s sequential -o <doc_directory>

It reads the GEOS TimeHistory output (sealedHeatingHistory.hdf5) of both runs and writes a single lightweight CSV,
the file kept under version control:

    sealedHeating.csv    temperature, exact and GEOS pressure, porosity and stresses of both couplings

It exits with a non-zero status if a relative difference exceeds TOLERANCE.

Exact solution
--------------
The sample is laterally confined (eps_xx = eps_yy = 0), free on top (total sigma_zz = 0) and sealed, without
gravity, and its temperature is uniform. With I( T ) the integral of alpha from the initial temperature to T, the
vertical strain, the Biot porosity and the fluid mass per unit volume are
    eps_zz = ( 3 K I + b p ) / ( K + 4 G / 3 )
    phi    = phi_0 + b eps_zz + ( b - phi_0 ) / K_s p - 3 ( b - phi_0 ) I
    phi rho_0 exp( c_f p - beta_f T ) = phi_0 rho_0          (no fluid leaves the sample)
The last equation is one scalar equation for p at each temperature, solved by Newton's method to round-off.
The effective stresses are sigma_xx = ( K - 2 G / 3 ) eps_zz - 3 K I and sigma_zz = ( K + 4 G / 3 ) eps_zz - 3 K I.

The column "pressure_constantPorosityCoefficient" is the exact pressure when the porosity uses a constant
coefficient alpha( 0 ), as it did before the porosity update used the SolidThermalExpansion model; it shows how far
the test is from that error.
"""
import argparse
import os
import sys

import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))

# --- case properties (ThermoPoroElasticSealedHeating_base.xml) ---
BULK_MODULUS = 9.12345679e9
SHEAR_MODULUS = 6.00813e9
GRAIN_BULK_MODULUS = 1.0e11
REFERENCE_POROSITY = 0.2
FLUID_COMPRESSIBILITY = 5e-10
FLUID_THERMAL_EXPANSION = 2e-4
T_INITIAL = 0.0
TABLE_KNOTS = np.array([0.0, 20.0, 50.0, 70.0])
TABLE_VALUES = np.array([0.2e-5, 0.5e-5, 0.6e-5, 1.6e-5])

TOLERANCE = {"fim": 1e-10, "sequential": 1e-8}   # the sequential runs converge to newtonTol only


def thermalStrain(T):
    """Integral of alpha from T_INITIAL to T, exact for the piecewise-linear table."""
    lo, hi = sorted((T_INITIAL, T))
    x = np.concatenate(([lo], TABLE_KNOTS[(TABLE_KNOTS > lo) & (TABLE_KNOTS < hi)], [hi]))
    alpha = np.interp(x, TABLE_KNOTS, TABLE_VALUES)
    value = np.sum(0.5 * (alpha[1:] + alpha[:-1]) * np.diff(x))
    return value if T >= T_INITIAL else -value


def exactState(T, porosityStrain=None):
    """Pressure, porosity and effective stresses sigma_xx, sigma_zz at temperature T.

    porosityStrain, if given, replaces the thermal strain in the porosity (used for the pre-fix comparison).
    """
    K, G, Ks, phi0 = BULK_MODULUS, SHEAR_MODULUS, GRAIN_BULK_MODULUS, REFERENCE_POROSITY
    b = 1.0 - K / Ks
    M = K + 4.0 * G / 3.0
    I = thermalStrain(T)
    Iphi = I if porosityStrain is None else porosityStrain

    def porosity(p):
        return phi0 + b * (3.0 * K * I + b * p) / M + (b - phi0) / Ks * p - 3.0 * (b - phi0) * Iphi

    p = 0.0
    for _ in range(100):
        residual = np.log(porosity(p) / phi0) + FLUID_COMPRESSIBILITY * p - FLUID_THERMAL_EXPANSION * (T - T_INITIAL)
        derivative = (b * b / M + (b - phi0) / Ks) / porosity(p) + FLUID_COMPRESSIBILITY
        step = residual / derivative
        p -= step
        if abs(step) <= 1e-15 * max(abs(p), 1.0):
            break
    ezz = (3.0 * K * I + b * p) / M
    return p, porosity(p), (K - 2.0 * G / 3.0) * ezz - 3.0 * K * I, M * ezz - 3.0 * K * I


def read(directory):
    with h5py.File(os.path.join(directory, "sealedHeatingHistory.hdf5"), "r") as h:
        temperature = np.array(h["temperature"])[:, 0]
        pressure = np.array(h["pressure"])
        porosity = np.array(h["rockPorosity_porosity"]).reshape(len(temperature), -1)
        stress = np.array(h["rockSolid_stress"])
    return temperature, pressure, porosity, stress[:, :, 0], stress[:, :, 2]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--fim", required=True, help="run directory of the fully implicit deck")
    parser.add_argument("-s", "--sequential", required=True, help="run directory of the sequential deck")
    parser.add_argument("-o", "--output", default=HERE, help="directory where the CSV is written")
    args = parser.parse_args()

    runs = {"fim": read(args.fim), "sequential": read(args.sequential)}
    temperature = runs["fim"][0]
    exact = np.array([exactState(T) for T in temperature])
    preFix = np.array([exactState(T, porosityStrain=TABLE_VALUES[0] * (T - T_INITIAL))[0] for T in temperature])

    ok = True
    for name, (T, p, phi, sxx, szz) in runs.items():
        assert np.allclose(T, temperature)
        errors = {"pressure": np.max(np.abs(p - exact[:, 0, None])) / np.max(np.abs(exact[:, 0])),
                  "porosity change": np.max(np.abs(phi - exact[:, 1, None])) / np.max(np.abs(exact[:, 1] - REFERENCE_POROSITY)),
                  "sigma_xx": np.max(np.abs(sxx - exact[:, 2, None])) / np.max(np.abs(exact[:, 2])),
                  "sigma_zz": np.max(np.abs(szz - exact[:, 3, None])) / np.max(np.abs(exact[:, 3]))}
        print(f"{name:10s} max relative difference GEOS / exact: "
              + ", ".join(f"{k} {v:.1e}" for k, v in errors.items()))
        ok = ok and max(errors.values()) <= TOLERANCE[name]
    print(f"a constant porosity coefficient alpha( 0 ) would change the final pressure by "
          f"{100 * (preFix[-1] - exact[-1, 0]) / exact[-1, 0]:.1f}%")

    header = ("temperature,pressure_analytical,porosity_analytical,sigma_xx_analytical,sigma_zz_analytical,"
              "pressure_constantPorosityCoefficient,"
              "pressure_fim,porosity_fim,sigma_xx_fim,pressure_sequential,porosity_sequential,sigma_xx_sequential")
    columns = [temperature, exact[:, 0], exact[:, 1], exact[:, 2], exact[:, 3], preFix]
    for name in ("fim", "sequential"):
        _, p, phi, sxx, _ = runs[name]
        columns += [p[:, 0], phi[:, 0], sxx[:, 0]]
    np.savetxt(os.path.join(args.output, "sealedHeating.csv"), np.column_stack(columns), delimiter=",",
               header=header, comments="", fmt="%.12e")
    print("PASS" if ok else "FAIL: a difference exceeds its tolerance")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
