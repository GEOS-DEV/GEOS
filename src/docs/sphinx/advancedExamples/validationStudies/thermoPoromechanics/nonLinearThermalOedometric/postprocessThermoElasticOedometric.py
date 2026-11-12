"""LOCAL post-processing of the oedometric tests with a temperature-dependent thermal expansion coefficient.

Run by hand, after running the three benchmarks, each in its own directory:

    geos -i ThermoElasticOedometric_linearTEC_benchmark.xml          # in a linear/ directory
    geos -i ThermoElasticOedometric_tecTable_benchmark.xml           # in a tabulated/ directory
    geos -i ThermoElasticOedometric_linearTEC_cooling_benchmark.xml  # in a cooling/ directory
    python3 postprocessThermoElasticOedometric.py -l linear -t tabulated -c cooling -o <doc_directory>

It reads the GEOS TimeHistory outputs (.hdf5) and writes a single lightweight CSV:

    thermoElasticOedometric.csv    temperature, GEOS stress and analytical stress of the three cases

That CSV is the file kept under version control. The .hdf5 files are not: the documentation
figure is produced by plotThermoElasticOedometric.py, which reads this CSV only.

Analytical solution
-------------------
The sample is fixed laterally (x and y) and free on its top face, without gravity or pore pressure,
hence a state of uniaxial strain:
    eps_xx = eps_yy = 0,  sigma_zz = 0
The thermo-elastic law, with the thermal strain eps_th = int_{T0}^{T} alpha(s) ds, then gives
    sigma_xx = sigma_yy = - E / ( 1 - nu ) * int_{T0}^{T} alpha(s) ds
where alpha(T) is linear (referenceDrainedLinearTEC, dDrainedLinearTEC_dTemperature, referenceTemperature) or piecewise linear
(TableFunction), and T0 is the initial temperature.
"""
import argparse
import os
import xml.etree.ElementTree as ElementTree

import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(HERE, "../../../../../../../inputFiles/thermoPoromechanics")


def elasticModuli(xmlFile):
    """Young's modulus and Poisson's ratio of the ElasticIsotropic model of a deck."""
    elastic = ElementTree.parse(xmlFile).find("Constitutive/ElasticIsotropic")
    K = float(elastic.get("defaultBulkModulus"))
    G = float(elastic.get("defaultShearModulus"))
    nu = (3.0 * K - 2.0 * G) / (6.0 * K + 2.0 * G)
    return 3.0 * K * (1.0 - 2.0 * nu), nu


def linearAlpha(xmlFile):
    """alpha(T) of a linear SolidThermalExpansion model, and its kinks (none)."""
    tec = ElementTree.parse(xmlFile).find("Constitutive/SolidThermalExpansion")
    a0 = float(tec.get("referenceDrainedLinearTEC"))
    slope = float(tec.get("dDrainedLinearTEC_dTemperature", "0"))
    T0 = float(tec.get("referenceTemperature", "0"))
    return (lambda T: a0 + slope * (np.asarray(T) - T0)), np.array([])


def tabulatedAlpha(xmlFile):
    """alpha(T) of a tabulated SolidThermalExpansion model (linear interpolation, constant outside), and its kinks."""
    root = ElementTree.parse(xmlFile)
    tableName = root.find("Constitutive/SolidThermalExpansion").get("drainedLinearTECTableName")
    table = [t for t in root.findall("Functions/TableFunction") if t.get("name") == tableName][0]
    parse = lambda s: np.array([float(v) for v in s.strip("{} ").split(",")])
    knots, values = parse(table.get("coordinates")), parse(table.get("values"))
    return (lambda T: np.interp(T, knots, values)), knots


def integral(alpha, kinks, a, b, n=201):
    """int_a^b alpha(s) ds for a piecewise linear alpha.

    The trapezoidal rule is exact on each linear piece, so it is exact up to round-off once the kinks of
    alpha lying in [a, b] are added to the grid.
    """
    lower, upper = min(a, b), max(a, b)
    s = np.union1d(np.linspace(lower, upper, n), kinks[(kinks > lower) & (kinks < upper)])
    if b < a:
        s = s[::-1]
    y = alpha(s)
    return np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(s))


def readCase(directory, stressFile, temperatureFile, stressKey):
    with h5py.File(os.path.join(directory, stressFile), "r") as h:
        stress = np.array(h[stressKey])
        time = np.array(h[stressKey + " Time"])[:, 0]
    with h5py.File(os.path.join(directory, temperatureFile), "r") as h:
        temperature = np.array(h["temperature"])[:, 0]
    return time, temperature, stress[:, 0, 0]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-l", "--linear", required=True, help="run directory of the linear heating benchmark")
    parser.add_argument("-t", "--tabulated", required=True, help="run directory of the tabulated heating benchmark")
    parser.add_argument("-c", "--cooling", required=True, help="run directory of the linear cooling benchmark")
    parser.add_argument("-o", "--output", default=HERE, help="directory where the CSV is written")
    args = parser.parse_args()

    linearDeck = os.path.join(INPUT_DIR, "ThermoElasticOedometric_linearTEC_base.xml")
    tabulatedDeck = os.path.join(INPUT_DIR, "ThermoElasticOedometric_tecTable_base.xml")
    coolingDeck = os.path.join(INPUT_DIR, "ThermoElasticOedometric_linearTEC_cooling_benchmark.xml")

    cases = [
        ("linear", args.linear, "stressHistory_rock.hdf5", "temperatureHistory_rock.hdf5",
         "rockSolid_linearTEC_stress", linearDeck, linearAlpha(linearDeck)),
        ("tabulated", args.tabulated, "tabulatedStressHistory_rock.hdf5", "tabulatedTemperatureHistory_rock.hdf5",
         "rockSolid_tecTable_stress", tabulatedDeck, tabulatedAlpha(tabulatedDeck)),
        ("cooling", args.cooling, "coolingStressHistory_rock.hdf5", "coolingTemperatureHistory_rock.hdf5",
         "rockSolid_linearTEC_stress", coolingDeck, linearAlpha(coolingDeck)),
    ]

    columns, header = [], []
    for name, directory, stressFile, temperatureFile, stressKey, deck, (alpha, kinks) in cases:
        time, temperature, sigma = readCase(directory, stressFile, temperatureFile, stressKey)
        E, nu = elasticModuli(deck)
        T0 = temperature[0]
        analytical = np.array([-E / (1.0 - nu) * integral(alpha, kinks, T0, T) for T in temperature])
        mask = np.abs(analytical) > 0.0
        print(f"{name:10s}: max relative difference GEOS / analytical = "
              f"{np.max(np.abs(sigma[mask] - analytical[mask]) / np.abs(analytical[mask])):.2e}")
        columns += [time, temperature, sigma, analytical]
        header += [f"time_{name}", f"temperature_{name}", f"sigma_xx_geos_{name}", f"sigma_xx_analytical_{name}"]

    output = os.path.join(args.output, "thermoElasticOedometric.csv")
    np.savetxt(output, np.column_stack(columns), delimiter=",", header=",".join(header), comments="", fmt="%.8e")
    print(f"written {output}")


if __name__ == "__main__":
    main()
