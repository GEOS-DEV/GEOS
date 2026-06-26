#!/usr/bin/env python3
"""Run and check a material-point VonMisesJ uniaxial-stress test.

The case prescribes strain increments in xx and uses stress control to keep
all other stress components at zero.  The compiled material-point driver output
is compared with the one-dimensional elastic-perfectly-plastic Von Mises
solution for the PFW ``verificationVonMises`` material.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import sys
from typing import Any, Dict, Iterable, List, Mapping, Sequence

import numpy as np

_THIS_DIR = Path(__file__).resolve().parents[2]
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import pfw_materials as material_db  # noqa: E402
from pfw_material_point import run_material_point  # noqa: E402

VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")


def _material_constants() -> Dict[str, float]:
    mat = material_db.verificationVonMises
    return {
        "E": float(mat["defaultYoungModulus"]),
        "nu": float(mat["defaultPoissonRatio"]),
        "K": float(mat["defaultBulkModulus"]),
        "G": float(mat["defaultShearModulus"]),
        "Y": float(mat["defaultYieldStrength"]),
        "rho0": float(mat["defaultDensity"]),
    }


def build_case(final_strain: float, n_steps: int, dt: float, stress_tolerance: float, fd_epsilon: float, max_iterations: int) -> Dict[str, Any]:
    """Return a compiled-driver case dictionary for uniaxial stress."""
    if n_steps <= 0:
        raise ValueError("n_steps must be positive")
    strain_increment = float(final_strain) / float(n_steps)
    constants = _material_constants()
    return {
        "name": "vonMisesJ_uniaxial_stress_material_point",
        "units": "dimensionless",
        "material": {"source": "pfw_materials.py", "name": "verificationVonMises"},
        "initial": {
            "stress": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "temperature": 300.0,
            "temperatureRate": 0.0,
            "specificInternalEnergy": 0.0,
            "density": constants["rho0"],
            "lengthScale": 1.0,
            "strengthScale": 1.0,
            "plasticStrain": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "materialDirection": {
            "type": "matrix",
            "rows": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "update": "fixed",
        },
        "time": {"dt": float(dt), "nSteps": int(n_steps)},
        "temperature": {"mode": "isothermal", "initial": 300.0},
        "energy": {"mode": "stressPower", "retentionFactor": 0.0, "outputAccumulatedStressPower": True},
        "solver": {
            "stressTolerance": float(stress_tolerance),
            "finiteDifferenceStrain": float(fd_epsilon),
            "maxIterations": int(max_iterations),
        },
        "control": [
            {"component": "xx", "mode": "strain", "value": strain_increment},
            {"component": "yy", "mode": "stress", "target": 0.0},
            {"component": "zz", "mode": "stress", "target": 0.0},
            {"component": "yz", "mode": "stress", "target": 0.0},
            {"component": "xz", "mode": "stress", "target": 0.0},
            {"component": "xy", "mode": "stress", "target": 0.0},
        ],
    }


def read_driver_csv(path: Path) -> List[Dict[str, float]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows: List[Dict[str, float]] = []
        for row in reader:
            parsed: Dict[str, float] = {}
            for key, value in row.items():
                if key is None:
                    continue
                if value is None or value == "":
                    parsed[key] = float("nan")
                else:
                    parsed[key] = float(value)
            rows.append(parsed)
    return rows


def expected_history(rows: Sequence[Mapping[str, float]], constants: Mapping[str, float]) -> Dict[str, np.ndarray]:
    """Return small-strain uniaxial-stress expected histories.

    For uniaxial tension with zero lateral stress, the Von Mises equivalent
    stress equals |sigma_xx|.  The perfect-plastic expected response is

      sigma_xx = min(E eps_xx, Y),
      eps_p_xx = max(0, eps_xx - Y/E),
      eps_p_yy = eps_p_zz = -0.5 eps_p_xx.
    """
    eps_xx = np.asarray([float(row["eps_xx"]) for row in rows], dtype=float)
    E = float(constants["E"])
    Y = float(constants["Y"])
    elastic_limit = Y / E
    plastic_axial = np.maximum(0.0, eps_xx - elastic_limit)
    sigma_xx = np.minimum(E * np.maximum(eps_xx, 0.0), Y)

    stress = np.zeros((len(rows), 6), dtype=float)
    stress[:, 0] = sigma_xx
    plastic = np.zeros((len(rows), 6), dtype=float)
    plastic[:, 0] = plastic_axial
    plastic[:, 1] = -0.5 * plastic_axial
    plastic[:, 2] = -0.5 * plastic_axial
    return {"eps_xx": eps_xx, "stress": stress, "plasticStrain": plastic}


def measured_array(rows: Sequence[Mapping[str, float]], prefix: str) -> np.ndarray:
    data = np.zeros((len(rows), 6), dtype=float)
    missing: List[str] = []
    for c, name in enumerate(VOIGT):
        key = f"{prefix}_{name}"
        if key not in rows[0]:
            missing.append(key)
            continue
        data[:, c] = [float(row[key]) for row in rows]
    if missing:
        raise KeyError("driver output is missing expected columns: " + ", ".join(missing))
    return data


def measured_plastic_array(rows: Sequence[Mapping[str, float]]) -> np.ndarray:
    direct = [f"plasticStrain_{i}" for i in range(6)]
    if all(key in rows[0] for key in direct):
        return np.asarray([[float(row[key]) for key in direct] for row in rows], dtype=float)
    # Some material wrappers may expose a scalar plastic-strain magnitude.  That is
    # not enough for this component-level check, so keep the error explicit.
    available = sorted(key for key in rows[0].keys() if "plastic" in key.lower())
    raise KeyError("driver output does not contain plasticStrain_0..plasticStrain_5 columns; available plastic columns: " + ", ".join(available))


def von_mises_equivalent(stress: np.ndarray) -> np.ndarray:
    sx, sy, sz, syz, sxz, sxy = [stress[:, i] for i in range(6)]
    return np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3.0 * (syz ** 2 + sxz ** 2 + sxy ** 2))


def write_expected_csv(path: Path, rows: Sequence[Mapping[str, float]], expected: Mapping[str, np.ndarray]) -> None:
    with path.open("w", newline="") as stream:
        fieldnames = ["step", "eps_xx"]
        fieldnames.extend(f"expected_stress_{name}" for name in VOIGT)
        fieldnames.extend(f"expected_plasticStrain_{i}" for i in range(6))
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for i, row in enumerate(rows):
            out: Dict[str, float] = {"step": float(row.get("step", i + 1)), "eps_xx": float(expected["eps_xx"][i])}
            for c, name in enumerate(VOIGT):
                out[f"expected_stress_{name}"] = float(expected["stress"][i, c])
            for c in range(6):
                out[f"expected_plasticStrain_{c}"] = float(expected["plasticStrain"][i, c])
            writer.writerow(out)


def _ensure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def plot_components(x: np.ndarray, measured: np.ndarray, expected: np.ndarray, ylabel: str, title: str, path: Path, component_labels: Sequence[str]) -> None:
    plt = _ensure_matplotlib()
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    for c, label in enumerate(component_labels):
        ax.plot(x, expected[:, c], linestyle="--", linewidth=1.2, label=f"expected {label}")
        ax.plot(x, measured[:, c], linewidth=1.0, label=f"driver {label}")
    ax.set_xlabel(r"prescribed $\epsilon_{xx}$")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize="x-small", ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_equivalent(x: np.ndarray, stress: np.ndarray, expected_stress: np.ndarray, constants: Mapping[str, float], path: Path) -> None:
    plt = _ensure_matplotlib()
    q = von_mises_equivalent(stress)
    q_expected = von_mises_equivalent(expected_stress)
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(x, q_expected, linestyle="--", linewidth=1.5, label="expected von Mises stress")
    ax.plot(x, q, linewidth=1.2, label="driver von Mises stress")
    ax.axhline(float(constants["Y"]), linestyle=":", linewidth=1.0, label="yield strength")
    ax.set_xlabel(r"prescribed $\epsilon_{xx}$")
    ax.set_ylabel("von Mises stress")
    ax.set_title("VonMisesJ material-point uniaxial-stress check")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def check_results(rows: Sequence[Mapping[str, float]], stress: np.ndarray, plastic: np.ndarray, expected: Mapping[str, np.ndarray], stress_tol: float, plastic_tol: float, lateral_stress_tol: float) -> Dict[str, Any]:
    stress_error = stress - expected["stress"]
    plastic_error = plastic - expected["plasticStrain"]
    lateral_indices = [1, 2, 3, 4, 5]
    converged_values = [int(round(float(row.get("converged", 1.0)))) for row in rows]
    summary: Dict[str, Any] = {
        "numSteps": len(rows),
        "finalEpsXx": float(expected["eps_xx"][-1]),
        "maxAbsStressError": float(np.nanmax(np.abs(stress_error))),
        "rmsStressError": float(math.sqrt(np.nanmean(stress_error ** 2))),
        "maxAbsPlasticStrainError": float(np.nanmax(np.abs(plastic_error))),
        "rmsPlasticStrainError": float(math.sqrt(np.nanmean(plastic_error ** 2))),
        "maxAbsStressControlResidual": float(np.nanmax(np.abs(stress[:, lateral_indices]))),
        "allConverged": all(value == 1 for value in converged_values),
        "finalMeasuredStress": {name: float(stress[-1, i]) for i, name in enumerate(VOIGT)},
        "finalExpectedStress": {name: float(expected["stress"][-1, i]) for i, name in enumerate(VOIGT)},
        "finalMeasuredPlasticStrain": {str(i): float(plastic[-1, i]) for i in range(6)},
        "finalExpectedPlasticStrain": {str(i): float(expected["plasticStrain"][-1, i]) for i in range(6)},
        "stressTolerance": float(stress_tol),
        "plasticStrainTolerance": float(plastic_tol),
        "lateralStressTolerance": float(lateral_stress_tol),
    }
    failures: List[str] = []
    if not summary["allConverged"]:
        failures.append("at least one driver step did not converge")
    if summary["maxAbsStressError"] > stress_tol:
        failures.append(f"stress error {summary['maxAbsStressError']:.6e} exceeds tolerance {stress_tol:.6e}")
    if summary["maxAbsPlasticStrainError"] > plastic_tol:
        failures.append(f"plastic-strain error {summary['maxAbsPlasticStrainError']:.6e} exceeds tolerance {plastic_tol:.6e}")
    if summary["maxAbsStressControlResidual"] > lateral_stress_tol:
        failures.append(f"stress-control residual {summary['maxAbsStressControlResidual']:.6e} exceeds tolerance {lateral_stress_tol:.6e}")
    summary["passed"] = not failures
    summary["failures"] = failures
    return summary


def postprocess(output_csv: Path, work_dir: Path, stress_tol: float, plastic_tol: float, lateral_stress_tol: float, make_plots: bool = True, check: bool = True) -> Dict[str, Any]:
    rows = read_driver_csv(output_csv)
    if not rows:
        raise RuntimeError(f"driver output CSV is empty: {output_csv}")
    constants = _material_constants()
    expected = expected_history(rows, constants)
    stress = measured_array(rows, "stress")
    plastic = measured_plastic_array(rows)

    expected_csv = work_dir / "vonMisesJ_uniaxial_stress_material_point_expected.csv"
    write_expected_csv(expected_csv, rows, expected)

    if make_plots:
        plot_components(
            expected["eps_xx"],
            stress,
            expected["stress"],
            "Cauchy stress",
            "Stress components: driver vs. expected",
            work_dir / "vonMisesJ_uniaxial_stress_material_point_stress.png",
            VOIGT,
        )
        plot_components(
            expected["eps_xx"],
            plastic,
            expected["plasticStrain"],
            "plastic strain component",
            "Plastic strain components: driver vs. expected",
            work_dir / "vonMisesJ_uniaxial_stress_material_point_plastic_strain.png",
            [str(i) for i in range(6)],
        )
        plot_equivalent(
            expected["eps_xx"],
            stress,
            expected["stress"],
            constants,
            work_dir / "vonMisesJ_uniaxial_stress_material_point_von_mises.png",
        )

    summary = check_results(rows, stress, plastic, expected, stress_tol, plastic_tol, lateral_stress_tol)
    summary.update({
        "outputCsv": str(output_csv),
        "expectedCsv": str(expected_csv),
        "stressPlot": str(work_dir / "vonMisesJ_uniaxial_stress_material_point_stress.png"),
        "plasticStrainPlot": str(work_dir / "vonMisesJ_uniaxial_stress_material_point_plastic_strain.png"),
        "vonMisesPlot": str(work_dir / "vonMisesJ_uniaxial_stress_material_point_von_mises.png"),
        "material": "verificationVonMises",
        "materialConstants": constants,
    })
    summary_path = work_dir / "vonMisesJ_uniaxial_stress_material_point_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if check and not summary["passed"]:
        raise RuntimeError("VonMisesJ material-point verification failed: " + "; ".join(summary["failures"]))
    return summary


def run_command_for_display(command: Sequence[Any]) -> str:
    return " ".join(str(part) for part in command)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", default=os.environ.get("GEOS_MPM_MATERIAL_POINT_DRIVER", "geos_mpm_material_point_driver"))
    parser.add_argument("--work-dir", default="vonMisesJ_uniaxial_stress_material_point")
    parser.add_argument("--final-strain", type=float, default=0.06)
    parser.add_argument("--n-steps", type=int, default=120)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--driver-stress-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--fd-epsilon", type=float, default=1.0e-8)
    parser.add_argument("--max-newton-iterations", type=int, default=25)
    parser.add_argument("--stress-tolerance", type=float, default=2.0e-4)
    parser.add_argument("--plastic-tolerance", type=float, default=5.0e-4)
    parser.add_argument("--lateral-stress-tolerance", type=float, default=5.0e-8)
    parser.add_argument("--prepare-only", action="store_true", help="write sidecar files but do not run the compiled driver")
    parser.add_argument("--no-check", action="store_true", help="do not fail on analytical-comparison tolerance violations")
    parser.add_argument("--no-plot", action="store_true", help="skip PNG plot generation")
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    case = build_case(
        final_strain=args.final_strain,
        n_steps=args.n_steps,
        dt=args.dt,
        stress_tolerance=args.driver_stress_tolerance,
        fd_epsilon=args.fd_epsilon,
        max_iterations=args.max_newton_iterations,
    )
    if not args.prepare_only and shutil.which(args.driver) is None and not Path(args.driver).is_file():
        raise RuntimeError(
            f"material-point driver executable not found: {args.driver!r}. "
            "Build it with setupMPM --enable-material-point-driver, or pass --driver /path/to/geos_mpm_material_point_driver."
        )

    result = run_material_point(
        case,
        executable=args.driver,
        work_dir=work_dir,
        dry_run=args.prepare_only,
        check=True,
    )

    print("material XML:", result.material_xml)
    print("load path CSV:", result.load_path_csv)
    print("output CSV:", result.output_csv)
    print("command:", run_command_for_display(result.command))

    if args.prepare_only:
        print("prepared material-point verification files only; rerun without --prepare-only to execute")
        return 0

    if not Path(result.output_csv).is_file():
        raise RuntimeError(f"driver completed but output CSV was not found: {result.output_csv}")

    summary = postprocess(
        Path(result.output_csv),
        work_dir,
        stress_tol=args.stress_tolerance,
        plastic_tol=args.plastic_tolerance,
        lateral_stress_tol=args.lateral_stress_tolerance,
        make_plots=not args.no_plot,
        check=not args.no_check,
    )
    print("summary JSON:", work_dir / "vonMisesJ_uniaxial_stress_material_point_summary.json")
    print("stress plot:", work_dir / "vonMisesJ_uniaxial_stress_material_point_stress.png")
    print("plastic strain plot:", work_dir / "vonMisesJ_uniaxial_stress_material_point_plastic_strain.png")
    print("von Mises plot:", work_dir / "vonMisesJ_uniaxial_stress_material_point_von_mises.png")
    print("passed:", summary["passed"])
    if summary.get("failures"):
        for failure in summary["failures"]:
            print("failure:", failure)
    return 0 if summary["passed"] or args.no_check else 1


if __name__ == "__main__":
    raise SystemExit(main())
