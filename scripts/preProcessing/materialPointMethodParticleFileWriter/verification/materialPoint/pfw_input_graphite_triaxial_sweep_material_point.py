#!/usr/bin/env python3
"""Run a Triaxial-style Graphite material-point pressure/orientation sweep.

The default sweep is a plane-strain, triaxial compression offive selected material
orientations, theta = 0, 30, 45, 60, and 90 degrees, each loaded at seven
background pressures, P0 = 0, 1, 3, 5, 10, 20, and 30 GPa.

GEOS uses tension-positive Cauchy stress, so a positive initial pressure P0 is
written as sigma_xx = sigma_yy = sigma_zz = -P0.  By default, the x direction
is held in plane strain, the y direction is stress controlled to -P0, and the
z direction is compressed at a prescribed true strain rate.  Passing
--transverse-control full-stress stress-controls both x and y to -P0, matching
the two independent transverse barostats used in the MD protocol more closely.
The plotted ordinate is the compression-positive axial stress convention,

    -sigma_zz - P0,

against the accumulated axial true strain epsilon_zz.

The case labelled P0 = 0 is treated as the paper's ambient-pressure case by
default: it is plotted and tagged as 0 GPa, but the actual stress-control target
is 1 atm unless --zero-pressure-control exact is requested.  In full transverse
stress-control mode, the default stress-control strain lower bound is xx=0,yy=0;
this is equivalent to enforcing a minimum lateral true strain of zero, or a
minimum lateral stretch of one, for the triaxial-compression analogue.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np

_THIS_DIR = Path(__file__).resolve().parents[2]
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from pfw_material_point import default_compiled_driver_executable, run_material_point  # noqa: E402
import pfw_materials  # noqa: E402

DEFAULT_THETA_DEG = (0.0, 30.0, 45.0, 60.0, 90.0)
DEFAULT_PRESSURE_GPA = (0.0, 1.0, 3.0, 5.0, 10.0, 20.0, 30.0)
VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")
DEFAULT_MATERIAL_NAME = "graphiteTriaxialTest"
DEFAULT_AMBIENT_PRESSURE_GPA = 1.01325e-4


def parse_float_list(text: str) -> List[float]:
    values = [float(item) for item in text.replace(";", ",").split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one numeric value")
    return values


def format_float_tag(value: float) -> str:
    """Return a filesystem-friendly, stable tag for a float."""
    if abs(value - round(value)) < 1.0e-12:
        return f"{int(round(value)):03d}"
    text = f"{value:.12g}".replace("-", "m").replace(".", "p")
    return text


def control_pressure_for_case(
    label_pressure_gpa: float,
    zero_pressure_control: str,
    ambient_pressure_gpa: float,
    regularized_zero_pressure_gpa: float,
    transverse_control: str = "plane-strain",
) -> float:
    """Return the pressure used in initial stresses and stress targets.

    The ambient-pressure case is denoted as 0 GPa, but the is 1 atm.  
    Keeping the plot label at 0 GPa while using a small positive stress-control 
    target avoids the singular exact-zero barostat target without changing the 
    plotted pressure scale in a visible way (not sure if the is important)

    ``auto`` keeps the physically closest ambient target for the default
    plane-strain-style protocol, but uses a slightly larger regularized target
    for full transverse stress control. We also limit strain to be nonnegative
    in transverse directions.
    """
    pressure = float(label_pressure_gpa)
    if abs(pressure) > 1.0e-14:
        return pressure
    mode = str(zero_pressure_control).strip().lower()
    if mode == "auto":
        mode = "regularized" if normalize_transverse_control(transverse_control) == "full-stress" else "ambient"
    if mode == "exact":
        return 0.0
    if mode == "ambient":
        return max(0.0, float(ambient_pressure_gpa))
    if mode == "regularized":
        return max(0.0, float(regularized_zero_pressure_gpa))
    raise ValueError(f"unknown zero-pressure control mode {zero_pressure_control!r}")


def response_reference_pressure(rows: Sequence[Mapping[str, float]], pressure_gpa: float) -> float:
    """Return the pressure that should be subtracted from axial stress plots."""
    if rows:
        try:
            return float(rows[0].get("_control_pressure_gpa", pressure_gpa))
        except (TypeError, ValueError):
            return float(pressure_gpa)
    return float(pressure_gpa)


def default_density_for_material(material_name: str) -> float:
    """Return the PFW default density for a material name, with a useful error."""
    material = getattr(pfw_materials, str(material_name), None)
    if not isinstance(material, Mapping):
        raise RuntimeError(f"could not find material {material_name!r} in pfw_materials.py")
    if "defaultDensity" not in material:
        raise RuntimeError(f"material {material_name!r} does not define defaultDensity")
    return float(material["defaultDensity"])


def pressure_elasticity_warning(material_name: str) -> Optional[str]:
    """Return a warning when the selected Graphite material has no elastic pressure dependence."""
    material = getattr(pfw_materials, str(material_name), None)
    if not isinstance(material, Mapping):
        return None
    if str(material.get("model", "")).lower() != "graphite":
        return None
    keys = (
        "defaultYoungModulusTransversePressureDerivative",
        "defaultYoungModulusAxialPressureDerivative",
        "defaultShearModulusAxialTransversePressureDerivative",
    )
    derivatives = [abs(float(material.get(key, 0.0))) for key in keys]
    if max(derivatives) <= 0.0:
        return (
            f"selected material {material_name!r} has zero pressure-derivative elastic moduli; "
            "curves at different confining pressures will mostly differ only by the subtracted initial pressure. "
            f"Use --material-name {DEFAULT_MATERIAL_NAME} for the pressure-sensitive Graphite validation material."
        )
    return None


def graphite_material_frame(theta_deg: float) -> List[List[float]]:
    """Return the row-wise GEOS material frame for the Triaxial theta convention.

    Row 0 is the Graphite basal-plane normal.  The lab z-axis is the compression
    direction.  theta = 0 deg puts compression within the layers, while
    theta = 90 deg makes compression normal to the layers.
    """
    theta = math.radians(float(theta_deg))
    normal = np.asarray([0.0, math.cos(theta), math.sin(theta)], dtype=float)
    tangent0 = np.asarray([1.0, 0.0, 0.0], dtype=float)
    tangent0 -= float(np.dot(tangent0, normal)) * normal
    tangent0 /= float(np.linalg.norm(tangent0))
    tangent1 = np.cross(normal, tangent0)
    tangent1 /= float(np.linalg.norm(tangent1))
    return [normal.tolist(), tangent0.tolist(), tangent1.tolist()]


TRANSVERSE_CONTROL_ALIASES = {
    "plane": "plane-strain",
    "plane-strain": "plane-strain",
    "plane_strain": "plane-strain",
    "planestrain": "plane-strain",
    "planestrainxstressy": "plane-strain",
    "strain-x-stress-y": "plane-strain",
    "strain_x_stress_y": "plane-strain",
    "one-stress": "plane-strain",
    "single-stress": "plane-strain",
    "full": "full-stress",
    "full-stress": "full-stress",
    "full_stress": "full-stress",
    "fullstress": "full-stress",
    "stressxstressy": "full-stress",
    "stress": "full-stress",
    "stress-stress": "full-stress",
    "stress_stress": "full-stress",
    "barostat": "full-stress",
    "barostats": "full-stress",
}


def normalize_transverse_control(value: str) -> str:
    try:
        return TRANSVERSE_CONTROL_ALIASES[str(value).strip().lower()]
    except KeyError as exc:
        valid = ", ".join(sorted({"plane-strain", "full-stress"}))
        raise argparse.ArgumentTypeError(f"invalid transverse control {value!r}; expected one of {valid}") from exc


def stress_controlled_components(transverse_control: str) -> List[str]:
    mode = normalize_transverse_control(transverse_control)
    if mode == "full-stress":
        return ["xx", "yy"]
    return ["yy"]


def normalize_optional_bound(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if text.lower() in ("", "none", "off"):
        return None
    return text


def default_stress_control_min_strain(transverse_control: str, user_value: Optional[str]) -> Optional[Any]:
    if user_value is not None:
        return normalize_optional_bound(user_value)
    if normalize_transverse_control(transverse_control) == "full-stress":
        return {"xx": 0.0, "yy": 0.0}
    return None


def bound_for_metadata(value: Optional[Any]) -> str:
    if value is None:
        return ""
    if isinstance(value, Mapping):
        return ",".join(f"{key}={val}" for key, val in value.items())
    return str(value)


def control_protocol_name(transverse_control: str) -> str:
    mode = normalize_transverse_control(transverse_control)
    if mode == "full-stress":
        return "stressControlXY_trueStrainRateZ"
    return "planeStrainX_stressControlY_trueStrainRateZ"


def build_control_protocol(transverse_control: str, target_stress: float, strain_rate_us: float) -> List[Dict[str, Any]]:
    mode = normalize_transverse_control(transverse_control)
    if mode == "full-stress":
        xx_control = {"component": "xx", "mode": "stress", "target": float(target_stress)}
    else:
        # Plane-strain analogue of the quasi-2D Triaxial setup: the x direction
        # is treated as the out-of-plane direction and held fixed while y is the
        # expanding periodic direction controlled to the confining pressure.
        xx_control = {"component": "xx", "mode": "strain", "value": 0.0}
    return [
        xx_control,
        {"component": "yy", "mode": "stress", "target": float(target_stress)},
        {"component": "zz", "mode": "trueStrainRate", "rate": float(strain_rate_us)},
        # The local driver does not model the LAMMPS cell-tilt kinematics, so keep
        # shear deformation components fixed at zero for this protocol-level check.
        {"component": "yz", "mode": "strain", "value": 0.0},
        {"component": "xz", "mode": "strain", "value": 0.0},
        {"component": "xy", "mode": "strain", "value": 0.0},
    ]


def build_case(
    theta_deg: float,
    pressure_gpa: float,
    control_pressure_gpa: float,
    dt: float,
    n_steps: int,
    strain_rate_us: float,
    material_name: str,
    initial_density: float,
    length_scale: float,
    strength_scale: float,
    stress_tolerance: float,
    fd_epsilon: float,
    max_iterations: int,
    max_line_search_iterations: int,
    max_stress_bracket_iterations: int,
    max_stress_bisection_iterations: int,
    stress_bracket_initial_scale: float,
    stress_bracket_max_strain: float,
    stress_bracket_growth: float,
    stress_control_algorithm: str,
    stress_control_regularization: float,
    stress_control_max_strain_correction: float,
    stress_control_servo_compliance: float,
    stress_control_servo_relaxation: float,
    stress_control_servo_derivative_floor: float,
    stress_control_servo_iterations: int,
    stress_control_pattern_iterations: int,
    stress_control_pattern_initial_step: float,
    stress_control_pattern_min_step: float,
    stress_control_pattern_shrink: float,
    stress_control_pattern_growth: float,
    stress_failure_policy: str,
    stress_control_diagnostics: bool,
    stress_control_diagnostics_level: str,
    transverse_control: str,
    stress_control_min_strain: Optional[Any],
    stress_control_max_strain: Optional[Any],
) -> Dict[str, Any]:
    if n_steps <= 0:
        raise ValueError("n_steps must be positive")
    pressure = float(pressure_gpa)
    control_pressure = float(control_pressure_gpa)
    target_stress = -control_pressure
    theta_tag = format_float_tag(float(theta_deg))
    pressure_tag = format_float_tag(pressure)
    name = f"graphite_triaxial_theta_{theta_tag}_p_{pressure_tag}_material_point"
    return {
        "name": name,
        "units": "mm_us_mg_GPa_K",
        "metadata": {
            "labelPressureGPa": pressure,
            "controlPressureGPa": control_pressure,
            "stressControlMinStrain": bound_for_metadata(stress_control_min_strain),
            "stressControlMaxStrain": bound_for_metadata(stress_control_max_strain),
        },
        "material": {"source": "pfw_materials.py", "name": material_name},
        "initial": {
            "stress": [target_stress, target_stress, target_stress, 0.0, 0.0, 0.0],
            "temperature": 300.0,
            "temperatureRate": 0.0,
            "specificInternalEnergy": 0.0,
            "density": float(initial_density),
            "lengthScale": float(length_scale),
            "strengthScale": float(strength_scale),
            "damage": 0.0,
            "basalPlaneDamage": 0.0,
            "comminutionDamage": 0.0,
        },
        "materialDirection": {
            "type": "matrix",
            "rows": graphite_material_frame(float(theta_deg)),
            "update": "graphite",
        },
        "time": {"dt": float(dt), "nSteps": int(n_steps)},
        "temperature": {"mode": "isothermal", "initial": 300.0},
        "energy": {"mode": "stressPower", "retentionFactor": 0.0, "outputAccumulatedStressPower": True},
        "solver": {
            "stressTolerance": float(stress_tolerance),
            "finiteDifferenceStrain": float(fd_epsilon),
            "maxIterations": int(max_iterations),
            "maxLineSearchIterations": int(max_line_search_iterations),
            "maxStressBracketIterations": int(max_stress_bracket_iterations),
            "maxStressBisectionIterations": int(max_stress_bisection_iterations),
            "stressBracketInitialScale": float(stress_bracket_initial_scale),
            "stressBracketMaxStrain": float(stress_bracket_max_strain),
            "stressBracketGrowth": float(stress_bracket_growth),
            "stressControlAlgorithm": str(stress_control_algorithm),
            "stressControlRegularization": float(stress_control_regularization),
            "stressControlMaxStrainCorrection": float(stress_control_max_strain_correction),
            "stressControlServoCompliance": float(stress_control_servo_compliance),
            "stressControlServoRelaxation": float(stress_control_servo_relaxation),
            "stressControlServoDerivativeFloor": float(stress_control_servo_derivative_floor),
            "stressControlServoIterations": int(stress_control_servo_iterations),
            "stressControlPatternIterations": int(stress_control_pattern_iterations),
            "stressControlPatternInitialStep": float(stress_control_pattern_initial_step),
            "stressControlPatternMinStep": float(stress_control_pattern_min_step),
            "stressControlPatternShrink": float(stress_control_pattern_shrink),
            "stressControlPatternGrowth": float(stress_control_pattern_growth),
            "stressControlFailurePolicy": str(stress_failure_policy),
            "stressControlMinStrain": stress_control_min_strain,
            "stressControlMaxStrain": stress_control_max_strain,
            "stressControlDiagnostics": {
                "enabled": bool(stress_control_diagnostics),
                "level": str(stress_control_diagnostics_level),
            },
        },
        "control": build_control_protocol(str(transverse_control), target_stress, float(strain_rate_us)),
        "output": {"file": f"{name}.csv"},
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
    if not rows:
        raise RuntimeError(f"driver output CSV is empty: {path}")
    missing = [key for key in ("eps_zz", "stress_zz", "stress_xx", "stress_yy", "converged") if key not in rows[0]]
    if missing:
        raise KeyError(f"driver output {path} is missing columns: {', '.join(missing)}")
    return rows


def case_work_dir(root: Path, theta_deg: float, pressure_gpa: float) -> Path:
    return root / "cases" / f"theta_{format_float_tag(theta_deg)}_p_{format_float_tag(pressure_gpa)}"


def run_or_prepare_case(
    case: Mapping[str, Any],
    driver: str,
    work_dir: Path,
    dry_run: bool,
    reuse_existing: bool,
) -> Tuple[Path, List[str], Optional[int], Path]:
    expected_output = work_dir / f"{case['name']}.csv"
    expected_diagnostics = work_dir / f"{case['name']}.stress_control_diagnostics.csv"
    if reuse_existing and expected_output.is_file() and not dry_run:
        return expected_output, [driver, "<reused-existing-output>"], 0, expected_diagnostics
    # The sweep is a data-generation/robustness test.  Let the caller inspect any
    # partial CSV written before a stress-control failure instead of raising here.
    result = run_material_point(case, executable=driver, work_dir=work_dir, dry_run=dry_run, check=False)
    diagnostic_path = Path(result.stress_control_diagnostics_csv) if result.stress_control_diagnostics_csv else expected_diagnostics
    return Path(result.output_csv), list(result.command), result.returncode, diagnostic_path


def prepend_initial_row(rows: Sequence[Mapping[str, float]], pressure_gpa: float) -> List[Dict[str, float]]:
    reference_pressure = response_reference_pressure(rows, pressure_gpa)
    first = dict(rows[0])
    initial: Dict[str, float] = {key: 0.0 for key in first.keys()}
    initial["step"] = 0.0
    initial["time"] = 0.0
    initial["dt"] = 0.0
    for name in VOIGT:
        initial[f"eps_{name}"] = 0.0
        initial[f"stress_{name}"] = 0.0
    initial["stress_xx"] = -reference_pressure
    initial["stress_yy"] = -reference_pressure
    initial["stress_zz"] = -reference_pressure
    initial["_control_pressure_gpa"] = reference_pressure
    initial["converged"] = 1.0
    return [initial] + [dict(row) for row in rows]


def axial_response(rows: Sequence[Mapping[str, float]], pressure_gpa: float) -> Tuple[np.ndarray, np.ndarray]:
    data = prepend_initial_row(rows, pressure_gpa)
    reference_pressure = response_reference_pressure(data, pressure_gpa)
    ezz = np.asarray([float(row["eps_zz"]) for row in data], dtype=float)
    axial = np.asarray([-float(row["stress_zz"]) - reference_pressure for row in data], dtype=float)
    return ezz, axial


def write_processed_csv(path: Path, records: Sequence[Mapping[str, Any]]) -> None:
    fieldnames = [
        "theta_deg",
        "pressure_gpa",
        "control_pressure_gpa",
        "step",
        "time",
        "eps_zz",
        "minus_sigma_zz_minus_p0",
        "stress_xx",
        "stress_yy",
        "stress_zz",
        "eps_xx",
        "eps_yy",
        "stress_control_min_strain_xx",
        "stress_control_min_strain_yy",
        "stress_control_max_strain_xx",
        "stress_control_max_strain_yy",
        "confining_stress_residual_xx",
        "confining_stress_residual_yy",
        "converged",
        "newtonIterations",
        "stressResidualNorm",
        "density",
        "jacobian",
        "effectiveBulkModulus",
        "effectiveShearModulus",
        "damage",
        "basalPlaneDamage",
        "comminutionDamage",
        "caseOutputCsv",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({key: record.get(key, "") for key in fieldnames})




def write_case_summary_csv(path: Path, records: Sequence[Mapping[str, Any]]) -> None:
    if not records:
        return
    fieldnames = [
        "thetaDeg",
        "pressureGPa",
        "controlPressureGPa",
        "transverseControl",
        "stressControlledComponents",
        "numSteps",
        "allConverged",
        "numNonConvergedSteps",
        "firstNonConvergedStep",
        "lastNonConvergedStep",
        "maxStressResidualNorm",
        "finalStressResidualNorm",
        "maxAbsConfiningStressResidual",
        "maxAbsConfiningStressResidualXx",
        "maxAbsConfiningStressResidualYy",
        "maxAbsPlaneStrainEpsXx",
        "maxAbsTransverseStrainXx",
        "maxAbsTransverseStrainYy",
        "maxAbsTransverseStrain",
        "minEpsXx",
        "minEpsYy",
        "stressControlMinStrainXx",
        "stressControlMinStrainYy",
        "stressControlMaxStrainXx",
        "stressControlMaxStrainYy",
        "minJacobian",
        "maxDensity",
        "finalEpsZz",
        "targetFinalEpsZz",
        "finalEpsZzError",
        "finalMinusSigmaZzMinusP0",
        "partialCurve",
        "passed",
        "warnings",
        "failures",
        "outputCsv",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            row = {key: record.get(key, "") for key in fieldnames}
            for key in ("stressControlledComponents", "warnings", "failures"):
                value = row.get(key, "")
                if isinstance(value, (list, tuple)):
                    row[key] = ";".join(str(item) for item in value)
            writer.writerow(row)


def _maybe_float(value: Any) -> Any:
    if value is None or value == "":
        return ""
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def read_stress_control_diagnostics(
    path: Path,
    theta_deg: float,
    pressure_gpa: float,
    case_name: str,
) -> List[Dict[str, Any]]:
    if not path.is_file():
        return []
    rows: List[Dict[str, Any]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            parsed: Dict[str, Any] = {
                "theta_deg": float(theta_deg),
                "pressure_gpa": float(pressure_gpa),
                "caseName": case_name,
                "diagnosticCsv": str(path),
            }
            for key, value in row.items():
                if key is not None:
                    parsed[key] = _maybe_float(value)
            rows.append(parsed)
    return rows


def write_diagnostics_csv(path: Path, records: Sequence[Mapping[str, Any]]) -> None:
    if not records:
        return
    fieldnames: List[str] = ["theta_deg", "pressure_gpa", "caseName", "diagnosticCsv"]
    for record in records:
        for key in record.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({key: record.get(key, "") for key in fieldnames})



def write_diagnostics_summary_csv(path: Path, records: Sequence[Mapping[str, Any]]) -> None:
    """Write a compact per-case stress-control diagnostic summary."""
    if not records:
        return
    grouped: Dict[Tuple[float, float, str], List[Mapping[str, Any]]] = {}
    for record in records:
        key = (
            float(record.get("theta_deg", float("nan"))),
            float(record.get("pressure_gpa", float("nan"))),
            str(record.get("caseName", "")),
        )
        grouped.setdefault(key, []).append(record)

    fieldnames = [
        "theta_deg",
        "pressure_gpa",
        "caseName",
        "rows",
        "max_step",
        "final_stage",
        "final_residualNorm",
        "min_abs_residual_yy",
        "max_abs_residual_yy",
        "min_unknown_yy",
        "max_unknown_yy",
        "first_effectiveBulkModulus",
        "last_effectiveBulkModulus",
        "first_effectiveShearModulus",
        "last_effectiveShearModulus",
        "first_damage",
        "last_damage",
        "last_basalPlaneDamage",
        "last_comminutionDamage",
        "last_density",
        "last_jacobian",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for (theta, pressure, case_name), rows in sorted(grouped.items()):
            def f(row: Mapping[str, Any], key: str) -> float:
                try:
                    return float(row.get(key, float("nan")))
                except (TypeError, ValueError):
                    return float("nan")

            last = rows[-1]
            residuals = [abs(f(row, "residual_yy")) for row in rows if math.isfinite(f(row, "residual_yy"))]
            unknowns = [f(row, "unknown_yy") for row in rows if math.isfinite(f(row, "unknown_yy"))]
            writer.writerow({
                "theta_deg": theta,
                "pressure_gpa": pressure,
                "caseName": case_name,
                "rows": len(rows),
                "max_step": max((f(row, "step") for row in rows), default=float("nan")),
                "final_stage": last.get("stage", ""),
                "final_residualNorm": f(last, "residualNorm"),
                "min_abs_residual_yy": min(residuals) if residuals else "",
                "max_abs_residual_yy": max(residuals) if residuals else "",
                "min_unknown_yy": min(unknowns) if unknowns else "",
                "max_unknown_yy": max(unknowns) if unknowns else "",
                "first_effectiveBulkModulus": f(rows[0], "effectiveBulkModulus"),
                "last_effectiveBulkModulus": f(last, "effectiveBulkModulus"),
                "first_effectiveShearModulus": f(rows[0], "effectiveShearModulus"),
                "last_effectiveShearModulus": f(last, "effectiveShearModulus"),
                "first_damage": f(rows[0], "damage"),
                "last_damage": f(last, "damage"),
                "last_basalPlaneDamage": f(last, "basalPlaneDamage"),
                "last_comminutionDamage": f(last, "comminutionDamage"),
                "last_density": f(last, "density"),
                "last_jacobian": f(last, "jacobian"),
            })


def _ensure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def plot_one_theta(
    theta_deg: float,
    pressure_results: Mapping[float, Sequence[Mapping[str, float]]],
    output_path: Path,
    log_y: bool,
    y_floor: float,
) -> None:
    plt = _ensure_matplotlib()
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    has_positive = False
    for pressure in sorted(pressure_results.keys()):
        ezz, axial = axial_response(pressure_results[pressure], pressure)
        y = np.where(axial > y_floor, axial, np.nan) if log_y else axial
        has_positive = has_positive or bool(np.any(np.isfinite(y) & (y > 0.0)))
        ax.plot(ezz, y, linewidth=1.2, label=f"{pressure:g} GPa")
    ax.set_xlabel(r"axial true strain, $\epsilon_{zz}$")
    ax.set_ylabel(r"$-\sigma_{zz} - P_0$ (GPa)")
    ax.set_title(rf"Graphite material-point sweep, $\theta={theta_deg:g}^\circ$")
    if log_y and has_positive:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize="small", title=r"$P_0$")
    ax.set_xlim(min(ax.get_xlim()[0], -0.5), 0.0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_five_panel(
    grouped_results: Mapping[float, Mapping[float, Sequence[Mapping[str, float]]]],
    theta_order: Sequence[float],
    output_path: Path,
    log_y: bool,
    y_floor: float,
) -> None:
    plt = _ensure_matplotlib()
    n = len(theta_order)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 4.6), sharex=True, sharey=True)
    if n == 1:
        axes = [axes]
    has_positive = False
    for ax, theta in zip(axes, theta_order):
        pressure_results = grouped_results.get(theta, {})
        for pressure in sorted(pressure_results.keys()):
            ezz, axial = axial_response(pressure_results[pressure], pressure)
            y = np.where(axial > y_floor, axial, np.nan) if log_y else axial
            has_positive = has_positive or bool(np.any(np.isfinite(y) & (y > 0.0)))
            ax.plot(ezz, y, linewidth=1.0, label=f"{pressure:g}")
        ax.set_title(rf"$\theta={theta:g}^\circ$")
        ax.set_xlabel(r"$\epsilon_{zz}$")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.5, 0.0)
    if log_y and has_positive:
        for ax in axes:
            ax.set_yscale("log")
    axes[0].set_ylabel(r"$-\sigma_{zz} - P_0$ (GPa)")
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(7, len(handles)), title=r"$P_0$ (GPa)")
    fig.suptitle("Triaxial-style Graphite material-point response")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_stress_component_time_five_panel(
    grouped_results: Mapping[float, Mapping[float, Sequence[Mapping[str, float]]]],
    theta_order: Sequence[float],
    component: str,
    output_path: Path,
    residual: bool = False,
) -> None:
    """Plot a stress-controlled component, or its target residual, versus time."""
    plt = _ensure_matplotlib()
    n = len(theta_order)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 4.6), sharex=True, sharey=True)
    if n == 1:
        axes = [axes]
    stress_key = f"stress_{component}"
    for ax, theta in zip(axes, theta_order):
        pressure_results = grouped_results.get(theta, {})
        for pressure in sorted(pressure_results.keys()):
            reference_pressure = response_reference_pressure(pressure_results[pressure], pressure)
            data = prepend_initial_row(pressure_results[pressure], pressure)
            time = np.asarray([float(row.get("time", 0.0)) for row in data], dtype=float)
            stress = np.asarray([float(row.get(stress_key, float("nan"))) for row in data], dtype=float)
            if residual:
                y = stress + reference_pressure
            else:
                y = stress
            ax.plot(time, y, linewidth=1.0, label=f"{pressure:g}")
            if time.size > 0:
                reference = 0.0 if residual else -reference_pressure
                ax.plot([float(time[0]), float(time[-1])], [reference, reference], linewidth=0.6, linestyle="--", alpha=0.35)
        if residual:
            ax.axhline(0.0, linewidth=0.8, linestyle="--")
        ax.set_title(rf"$\theta={theta:g}^\circ$")
        ax.set_xlabel(r"time ($\mu$s)")
        ax.grid(True, alpha=0.3)
    if residual:
        axes[0].set_ylabel(rf"$\sigma_{{{component}}} + P_0$ (GPa)")
        title = rf"Stress-control residual for $\sigma_{{{component}}}$"
    else:
        axes[0].set_ylabel(rf"$\sigma_{{{component}}}$ (GPa)")
        title = rf"Stress-controlled component $\sigma_{{{component}}}$ versus time"
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(7, len(handles)), title=r"$P_0$ (GPa)")
    fig.suptitle(title)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_output_field_time_five_panel(
    grouped_results: Mapping[float, Mapping[float, Sequence[Mapping[str, float]]]],
    theta_order: Sequence[float],
    field: str,
    ylabel: str,
    output_path: Path,
) -> bool:
    """Plot an optional driver output field versus time and return whether data existed."""
    has_data = False
    plt = _ensure_matplotlib()
    n = len(theta_order)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 4.6), sharex=True, sharey=True)
    if n == 1:
        axes = [axes]
    for ax, theta in zip(axes, theta_order):
        pressure_results = grouped_results.get(theta, {})
        for pressure in sorted(pressure_results.keys()):
            data = prepend_initial_row(pressure_results[pressure], pressure)
            if field not in data[-1]:
                continue
            time = np.asarray([float(row.get("time", 0.0)) for row in data], dtype=float)
            values = np.asarray([float(row.get(field, float("nan"))) for row in data], dtype=float)
            good = np.isfinite(time) & np.isfinite(values)
            if np.any(good):
                has_data = True
                ax.plot(time[good], values[good], linewidth=1.0, label=f"{pressure:g}")
        ax.set_title(rf"$\theta={theta:g}^\circ$")
        ax.set_xlabel(r"time ($\mu$s)")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(ylabel)
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(7, len(handles)), title=r"$P_0$ (GPa)")
    fig.suptitle(f"Graphite material-point {field} versus time")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    if has_data:
        fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return has_data




def plot_stress_control_trial_residual_five_panel(
    diagnostic_rows: Sequence[Mapping[str, Any]],
    theta_order: Sequence[float],
    component: str,
    output_path: Path,
) -> None:
    """Plot stress-control residual sampled at every trial evaluation."""
    if not diagnostic_rows:
        return
    plt = _ensure_matplotlib()
    n = len(theta_order)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 4.6), sharex=False, sharey=True)
    if n == 1:
        axes = [axes]
    residual_key = f"residual_{component}"
    unknown_key = f"unknown_{component}"
    for ax, theta in zip(axes, theta_order):
        theta_rows = [row for row in diagnostic_rows if abs(float(row.get("theta_deg", float("nan"))) - float(theta)) < 1.0e-10]
        pressures = sorted({float(row.get("pressure_gpa", float("nan"))) for row in theta_rows if math.isfinite(float(row.get("pressure_gpa", float("nan"))))})
        for pressure in pressures:
            rows = [row for row in theta_rows if abs(float(row.get("pressure_gpa", float("nan"))) - pressure) < 1.0e-10]
            x = np.asarray([float(row.get(unknown_key, float("nan"))) for row in rows], dtype=float)
            y = np.asarray([float(row.get(residual_key, float("nan"))) for row in rows], dtype=float)
            good = np.isfinite(x) & np.isfinite(y)
            if np.any(good):
                ax.plot(x[good], y[good], linestyle="None", marker=".", markersize=2.5, alpha=0.55, label=f"{pressure:g}")
        ax.axhline(0.0, linewidth=0.8, linestyle="--")
        ax.set_title(rf"$\theta={theta:g}^\circ$")
        ax.set_xlabel(rf"trial $\Delta \epsilon_{{{component}}}$")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(rf"trial residual $\sigma_{{{component}}}-\sigma_{{{component}}}^{{target}}$ (GPa)")
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(7, len(handles)), title=r"$P_0$ (GPa)")
    fig.suptitle(rf"Stress-control trial residuals for $\sigma_{{{component}}}$")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def resolve_component_list(text: str, controlled_components: Sequence[str]) -> List[str]:
    normalized = str(text).strip().lower()
    if normalized in ("all", "*"):
        return list(VOIGT)
    if normalized in ("controlled", "stress-controlled", "stress_controlled"):
        return list(controlled_components)
    components = [item.strip().lower() for item in str(text).replace(";", ",").split(",") if item.strip()]
    invalid = [item for item in components if item not in VOIGT]
    if invalid:
        raise argparse.ArgumentTypeError(f"invalid stress component(s): {','.join(invalid)}; expected one of {','.join(VOIGT)}")
    return components


def summarize_case(
    rows: Sequence[Mapping[str, float]],
    theta_deg: float,
    pressure_gpa: float,
    control_pressure_gpa: float,
    output_csv: Path,
    target_final_ezz: float,
    stress_control_tolerance: float,
    final_strain_tolerance: float,
    plane_strain_tolerance: float,
    transverse_control: str,
    strict_convergence: bool = False,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    converged = [int(round(float(row.get("converged", 1.0)))) for row in rows]
    nonconverged_steps = [int(round(float(row.get("step", i + 1)))) for i, row in enumerate(rows) if converged[i] != 1]
    stress_residual_norms = np.asarray([float(row.get("stressResidualNorm", 0.0)) for row in rows], dtype=float)
    finite_stress_residual_norms = stress_residual_norms[np.isfinite(stress_residual_norms)]
    max_stress_residual_norm = float(np.nanmax(finite_stress_residual_norms)) if finite_stress_residual_norms.size else float("nan")
    final_stress_residual_norm = float(stress_residual_norms[-1]) if stress_residual_norms.size else float("nan")
    target_stress = -float(control_pressure_gpa)
    controlled = stress_controlled_components(transverse_control)
    stress_by_component = {name: np.asarray([float(row[f"stress_{name}"]) for row in rows], dtype=float) for name in controlled}
    eps_xx = np.asarray([float(row["eps_xx"]) for row in rows], dtype=float)
    eps_yy = np.asarray([float(row["eps_yy"]) for row in rows], dtype=float)
    ezz = np.asarray([float(row["eps_zz"]) for row in rows], dtype=float)

    def first_finite_value(field_name: str) -> float:
        for row in rows:
            try:
                value = float(row.get(field_name, float("nan")))
            except (TypeError, ValueError):
                value = float("nan")
            if math.isfinite(value):
                return value
        return float("nan")

    stress_min_xx = first_finite_value("stressControlStrainMin_xx")
    stress_min_yy = first_finite_value("stressControlStrainMin_yy")
    stress_max_xx = first_finite_value("stressControlStrainMax_xx")
    stress_max_yy = first_finite_value("stressControlStrainMax_yy")

    def finite_optional_array(field_name: str) -> np.ndarray:
        values: List[float] = []
        for row in rows:
            try:
                value = float(row.get(field_name, float("nan")))
            except (TypeError, ValueError):
                value = float("nan")
            if math.isfinite(value):
                values.append(value)
        return np.asarray(values, dtype=float)

    jacobian_values = finite_optional_array("jacobian")
    density_values = finite_optional_array("density")
    residual_by_component = {name: values - target_stress for name, values in stress_by_component.items()}
    max_confining_residual_by_component = {
        name: float(np.nanmax(np.abs(values))) for name, values in residual_by_component.items()
    }
    max_confining_residual = max(max_confining_residual_by_component.values(), default=0.0)
    max_plane_strain_error = float(np.nanmax(np.abs(eps_xx)))
    max_abs_transverse_strain_xx = float(np.nanmax(np.abs(eps_xx))) if eps_xx.size else float("nan")
    max_abs_transverse_strain_yy = float(np.nanmax(np.abs(eps_yy))) if eps_yy.size else float("nan")
    max_abs_transverse_strain = float(np.nanmax([max_abs_transverse_strain_xx, max_abs_transverse_strain_yy]))
    min_eps_xx = float(np.nanmin(eps_xx)) if eps_xx.size else float("nan")
    min_eps_yy = float(np.nanmin(eps_yy)) if eps_yy.size else float("nan")
    min_jacobian = float(np.nanmin(jacobian_values)) if jacobian_values.size else float("nan")
    max_density = float(np.nanmax(density_values)) if density_values.size else float("nan")
    final_ezz_error = float(abs(ezz[-1] - target_final_ezz))
    partial_curve = final_ezz_error > final_strain_tolerance
    failures: List[str] = []
    warnings: List[str] = []

    def convergence_issue(message: str) -> None:
        if strict_convergence:
            failures.append(message)
        else:
            warnings.append(message)

    if not all(value == 1 for value in converged):
        convergence_issue("at least one step did not converge")
    if max_confining_residual > stress_control_tolerance:
        convergence_issue(f"confining stress residual {max_confining_residual:.6e} exceeds tolerance {stress_control_tolerance:.6e}")
    if normalize_transverse_control(transverse_control) == "plane-strain" and max_plane_strain_error > plane_strain_tolerance:
        failures.append(f"plane-strain eps_xx error {max_plane_strain_error:.6e} exceeds tolerance {plane_strain_tolerance:.6e}")
    if normalize_transverse_control(transverse_control) == "full-stress":
        if math.isfinite(stress_min_xx) and math.isfinite(min_eps_xx) and min_eps_xx < stress_min_xx - 1.0e-10:
            failures.append(f"eps_xx minimum {min_eps_xx:.6e} fell below stress-control strain lower bound {stress_min_xx:.6e}")
        if math.isfinite(stress_min_yy) and math.isfinite(min_eps_yy) and min_eps_yy < stress_min_yy - 1.0e-10:
            failures.append(f"eps_yy minimum {min_eps_yy:.6e} fell below stress-control strain lower bound {stress_min_yy:.6e}")
        if math.isfinite(max_abs_transverse_strain) and max_abs_transverse_strain > 2.0:
            warnings.append(f"large transverse strain magnitude {max_abs_transverse_strain:.6e}; near-zero pressure target may be ill-conditioned")
        if math.isfinite(min_jacobian) and min_jacobian < 1.0e-6:
            warnings.append(f"minimum jacobian {min_jacobian:.6e}; stress-controlled cell is approaching a singular deformation")
    if partial_curve:
        failures.append(f"partial curve: final eps_zz error {final_ezz_error:.6e} exceeds tolerance {final_strain_tolerance:.6e}")

    case_summary: Dict[str, Any] = {
        "thetaDeg": float(theta_deg),
        "pressureGPa": float(pressure_gpa),
        "controlPressureGPa": float(control_pressure_gpa),
        "outputCsv": str(output_csv),
        "numSteps": len(rows),
        "allConverged": all(value == 1 for value in converged),
        "numNonConvergedSteps": len(nonconverged_steps),
        "firstNonConvergedStep": nonconverged_steps[0] if nonconverged_steps else None,
        "lastNonConvergedStep": nonconverged_steps[-1] if nonconverged_steps else None,
        "maxStressResidualNorm": max_stress_residual_norm,
        "finalStressResidualNorm": final_stress_residual_norm,
        "partialCurve": bool(partial_curve),
        "transverseControl": normalize_transverse_control(transverse_control),
        "stressControlledComponents": controlled,
        "maxAbsConfiningStressResidual": max_confining_residual,
        "maxAbsConfiningStressResidualByComponent": max_confining_residual_by_component,
        "maxAbsConfiningStressResidualYy": max_confining_residual_by_component.get("yy", float("nan")),
        "maxAbsConfiningStressResidualXx": max_confining_residual_by_component.get("xx", float("nan")),
        "maxAbsPlaneStrainEpsXx": max_plane_strain_error,
        "maxAbsTransverseStrainXx": max_abs_transverse_strain_xx,
        "maxAbsTransverseStrainYy": max_abs_transverse_strain_yy,
        "maxAbsTransverseStrain": max_abs_transverse_strain,
        "minEpsXx": min_eps_xx,
        "minEpsYy": min_eps_yy,
        "stressControlMinStrainXx": stress_min_xx,
        "stressControlMinStrainYy": stress_min_yy,
        "stressControlMaxStrainXx": stress_max_xx,
        "stressControlMaxStrainYy": stress_max_yy,
        "minJacobian": min_jacobian,
        "maxDensity": max_density,
        "finalEpsZz": float(ezz[-1]),
        "targetFinalEpsZz": float(target_final_ezz),
        "finalEpsZzError": final_ezz_error,
        "finalMinusSigmaZzMinusP0": float(-float(rows[-1]["stress_zz"]) - float(control_pressure_gpa)),
        "passed": not failures,
        "failures": failures,
        "warnings": warnings,
    }

    processed_rows: List[Dict[str, Any]] = []
    for row in prepend_initial_row(rows, pressure_gpa):
        processed_rows.append({
            "theta_deg": float(theta_deg),
            "pressure_gpa": float(pressure_gpa),
            "control_pressure_gpa": float(control_pressure_gpa),
            "step": int(round(float(row.get("step", 0.0)))),
            "time": float(row.get("time", 0.0)),
            "eps_zz": float(row.get("eps_zz", 0.0)),
            "minus_sigma_zz_minus_p0": -float(row.get("stress_zz", -control_pressure_gpa)) - float(control_pressure_gpa),
            "stress_xx": float(row.get("stress_xx", float("nan"))),
            "stress_yy": float(row.get("stress_yy", float("nan"))),
            "stress_zz": float(row.get("stress_zz", float("nan"))),
            "eps_xx": float(row.get("eps_xx", float("nan"))),
            "eps_yy": float(row.get("eps_yy", float("nan"))),
            "stress_control_min_strain_xx": float(row.get("stressControlStrainMin_xx", float("nan"))),
            "stress_control_min_strain_yy": float(row.get("stressControlStrainMin_yy", float("nan"))),
            "stress_control_max_strain_xx": float(row.get("stressControlStrainMax_xx", float("nan"))),
            "stress_control_max_strain_yy": float(row.get("stressControlStrainMax_yy", float("nan"))),
            "confining_stress_residual_xx": float(row.get("stress_xx", float("nan"))) - target_stress,
            "confining_stress_residual_yy": float(row.get("stress_yy", float("nan"))) - target_stress,
            "converged": int(round(float(row.get("converged", 1.0)))),
            "newtonIterations": int(round(float(row.get("newtonIterations", 0.0)))),
            "stressResidualNorm": float(row.get("stressResidualNorm", 0.0)),
            "density": row.get("density", ""),
            "jacobian": row.get("jacobian", ""),
            "effectiveBulkModulus": row.get("effectiveBulkModulus", ""),
            "effectiveShearModulus": row.get("effectiveShearModulus", ""),
            "damage": row.get("damage", ""),
            "basalPlaneDamage": row.get("basalPlaneDamage", ""),
            "comminutionDamage": row.get("comminutionDamage", ""),
            "caseOutputCsv": str(output_csv),
        })
    return case_summary, processed_rows


def run_sweep(args: argparse.Namespace) -> Dict[str, Any]:
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    theta_values = list(args.theta_deg)
    pressure_values = list(args.pressure_gpa)
    target_final_ezz = float(args.strain_rate) * float(args.dt) * int(args.n_steps)
    initial_density = float(args.initial_density) if args.initial_density is not None else default_density_for_material(args.material_name)
    args.transverse_control = normalize_transverse_control(args.transverse_control)
    controlled_components = stress_controlled_components(args.transverse_control)
    args.stress_diagnostic_components = resolve_component_list(args.stress_diagnostic_components, controlled_components)
    stress_control_min_strain = default_stress_control_min_strain(args.transverse_control, args.stress_control_min_strain)
    stress_control_max_strain = normalize_optional_bound(args.stress_control_max_strain)

    material_pressure_warning = pressure_elasticity_warning(args.material_name)
    if material_pressure_warning:
        print(f"warning: {material_pressure_warning}", file=sys.stderr)

    if not args.prepare_only and not args.postprocess_only:
        if shutil.which(args.driver) is None and not Path(args.driver).is_file():
            raise RuntimeError(
                f"material-point driver executable not found: {args.driver!r}. "
                "Build it with setupMPM --enable-material-point-driver, or pass --driver /path/to/geos_mpm_material_point_driver."
            )

    grouped_results: Dict[float, Dict[float, List[Dict[str, float]]]] = {}
    processed_rows: List[Dict[str, Any]] = []
    diagnostic_rows: List[Dict[str, Any]] = []
    case_summaries: List[Dict[str, Any]] = []
    prepared_cases: List[Dict[str, Any]] = []
    commands: List[Dict[str, Any]] = []
    failures: List[str] = []
    warnings: List[str] = []
    if material_pressure_warning:
        warnings.append(material_pressure_warning)

    for theta in theta_values:
        grouped_results.setdefault(float(theta), {})
        for pressure in pressure_values:
            control_pressure = control_pressure_for_case(
                float(pressure),
                args.zero_pressure_control,
                args.ambient_pressure_gpa,
                args.regularized_zero_pressure_gpa,
                args.transverse_control,
            )
            case = build_case(
                theta_deg=float(theta),
                pressure_gpa=float(pressure),
                control_pressure_gpa=control_pressure,
                dt=args.dt,
                n_steps=args.n_steps,
                strain_rate_us=args.strain_rate,
                material_name=args.material_name,
                initial_density=initial_density,
                length_scale=args.length_scale,
                strength_scale=args.strength_scale,
                stress_tolerance=args.driver_stress_tolerance,
                fd_epsilon=args.fd_epsilon,
                max_iterations=args.max_newton_iterations,
                max_line_search_iterations=args.max_line_search_iterations,
                max_stress_bracket_iterations=args.max_stress_bracket_iterations,
                max_stress_bisection_iterations=args.max_stress_bisection_iterations,
                stress_bracket_initial_scale=args.stress_bracket_initial_scale,
                stress_bracket_max_strain=args.stress_bracket_max_strain,
                stress_bracket_growth=args.stress_bracket_growth,
                stress_control_algorithm=args.driver_stress_control_algorithm,
                stress_control_regularization=args.stress_control_regularization,
                stress_control_max_strain_correction=args.stress_control_max_strain_correction,
                stress_control_servo_compliance=args.stress_control_servo_compliance,
                stress_control_servo_relaxation=args.stress_control_servo_relaxation,
                stress_control_servo_derivative_floor=args.stress_control_servo_derivative_floor,
                stress_control_servo_iterations=args.stress_control_servo_iterations,
                stress_control_pattern_iterations=args.stress_control_pattern_iterations,
                stress_control_pattern_initial_step=args.stress_control_pattern_initial_step,
                stress_control_pattern_min_step=args.stress_control_pattern_min_step,
                stress_control_pattern_shrink=args.stress_control_pattern_shrink,
                stress_control_pattern_growth=args.stress_control_pattern_growth,
                stress_failure_policy=args.driver_stress_failure_policy,
                stress_control_diagnostics=(not args.no_driver_stress_control_diagnostics),
                stress_control_diagnostics_level=args.driver_stress_control_diagnostics_level,
                transverse_control=args.transverse_control,
                stress_control_min_strain=stress_control_min_strain,
                stress_control_max_strain=stress_control_max_strain,
            )
            cdir = case_work_dir(work_dir, float(theta), float(pressure))
            prepared_cases.append({
                "thetaDeg": float(theta),
                "pressureGPa": float(pressure),
                "controlPressureGPa": float(control_pressure),
                "stressControlMinStrain": bound_for_metadata(stress_control_min_strain),
                "stressControlMaxStrain": bound_for_metadata(stress_control_max_strain),
                "workDir": str(cdir),
                "caseName": case["name"],
            })
            try:
                if args.postprocess_only:
                    output_csv = cdir / f"{case['name']}.csv"
                    diagnostic_csv = cdir / f"{case['name']}.stress_control_diagnostics.csv"
                    command = [args.driver, "<postprocess-only>"]
                    returncode = None
                else:
                    output_csv, command, returncode, diagnostic_csv = run_or_prepare_case(
                        case=case,
                        driver=args.driver,
                        work_dir=cdir,
                        dry_run=args.prepare_only,
                        reuse_existing=args.reuse_existing,
                    )
                commands.append({
                    "thetaDeg": float(theta),
                    "pressureGPa": float(pressure),
                    "controlPressureGPa": float(control_pressure),
                    "stressControlMinStrain": bound_for_metadata(stress_control_min_strain),
                    "stressControlMaxStrain": bound_for_metadata(stress_control_max_strain),
                    "outputCsv": str(output_csv),
                    "diagnosticCsv": str(diagnostic_csv),
                    "returncode": returncode,
                    "command": " ".join(str(part) for part in command),
                })
                if returncode not in (None, 0):
                    failures.append(f"theta={theta:g}, P={pressure:g}: driver returned {returncode}; attempting to use partial output")
                if args.prepare_only:
                    continue
                if not args.prepare_only:
                    diagnostic_rows.extend(read_stress_control_diagnostics(diagnostic_csv, float(theta), float(pressure), case["name"]))
                rows = read_driver_csv(output_csv)
                for row in rows:
                    row["_control_pressure_gpa"] = float(control_pressure)
                if not rows:
                    message = "driver produced no material-point rows; see the case .driver.log file"
                    failures.append(f"theta={theta:g}, P={pressure:g}: {message}")
                    case_summaries.append({
                        "thetaDeg": float(theta),
                        "pressureGPa": float(pressure),
                        "outputCsv": str(output_csv),
                        "numSteps": 0,
                        "partialCurve": True,
                        "passed": False,
                        "failures": [message],
                    })
                    continue
                grouped_results[float(theta)][float(pressure)] = rows
                case_summary, rows_for_csv = summarize_case(
                    rows,
                    theta_deg=float(theta),
                    pressure_gpa=float(pressure),
                    control_pressure_gpa=float(control_pressure),
                    output_csv=output_csv,
                    target_final_ezz=target_final_ezz,
                    stress_control_tolerance=args.stress_control_tolerance,
                    final_strain_tolerance=args.final_strain_tolerance,
                    plane_strain_tolerance=args.plane_strain_tolerance,
                    transverse_control=args.transverse_control,
                    strict_convergence=args.strict,
                )
                case_summaries.append(case_summary)
                processed_rows.extend(rows_for_csv)
                if case_summary.get("warnings"):
                    warnings.append(f"theta={theta:g}, P={pressure:g}: " + "; ".join(case_summary["warnings"]))
                if not case_summary["passed"]:
                    failures.append(f"theta={theta:g}, P={pressure:g}: " + "; ".join(case_summary["failures"]))
            except Exception as exc:  # keep the case label in the suite summary
                message = f"theta={theta:g}, P={pressure:g}: {exc}"
                failures.append(message)
                case_summaries.append({
                    "thetaDeg": float(theta),
                    "pressureGPa": float(pressure),
                    "passed": False,
                    "failures": [str(exc)],
                })
                if args.fail_fast:
                    raise

    plots: Dict[str, Any] = {}
    if not args.prepare_only and not args.no_plot and processed_rows:
        for theta in theta_values:
            plot_path = work_dir / f"graphite_triaxial_response_theta_{format_float_tag(float(theta))}.png"
            if grouped_results.get(float(theta)):
                plot_one_theta(float(theta), grouped_results[float(theta)], plot_path, log_y=not args.linear_y, y_floor=args.y_floor)
                plots[f"theta_{float(theta):g}"] = str(plot_path)
        theta_order = [float(theta) for theta in theta_values]
        combined = work_dir / "graphite_triaxial_response_5panel.png"
        plot_five_panel(grouped_results, theta_order, combined, log_y=not args.linear_y, y_floor=args.y_floor)
        plots["combined5Panel"] = str(combined)

        for component in args.stress_diagnostic_components:
            stress_plot = work_dir / f"graphite_triaxial_stress_{component}_vs_time_5panel.png"
            plot_stress_component_time_five_panel(grouped_results, theta_order, component, stress_plot, residual=False)
            plots[f"stress_{component}_vs_time_5panel"] = str(stress_plot)
            residual_plot = work_dir / f"graphite_triaxial_stress_residual_{component}_vs_time_5panel.png"
            plot_stress_component_time_five_panel(grouped_results, theta_order, component, residual_plot, residual=True)
            plots[f"stress_residual_{component}_vs_time_5panel"] = str(residual_plot)
            if diagnostic_rows:
                trial_plot = work_dir / f"graphite_triaxial_stress_control_trial_residual_{component}_5panel.png"
                plot_stress_control_trial_residual_five_panel(diagnostic_rows, theta_order, component, trial_plot)
                plots[f"stress_control_trial_residual_{component}_5panel"] = str(trial_plot)

        modulus_plots = (
            ("effectiveBulkModulus", r"effective bulk modulus (GPa)"),
            ("effectiveShearModulus", r"effective shear modulus (GPa)"),
        )
        for field, ylabel in modulus_plots:
            field_plot = work_dir / f"graphite_triaxial_{field}_vs_time_5panel.png"
            if plot_output_field_time_five_panel(grouped_results, theta_order, field, ylabel, field_plot):
                plots[f"{field}_vs_time_5panel"] = str(field_plot)

    processed_csv = work_dir / "graphite_triaxial_sweep_processed.csv"
    if processed_rows:
        write_processed_csv(processed_csv, processed_rows)

    diagnostics_csv = work_dir / "graphite_triaxial_stress_control_diagnostics.csv"
    diagnostics_summary_csv = work_dir / "graphite_triaxial_stress_control_diagnostics_summary.csv"
    if diagnostic_rows:
        write_diagnostics_csv(diagnostics_csv, diagnostic_rows)
        write_diagnostics_summary_csv(diagnostics_summary_csv, diagnostic_rows)

    case_summary_csv = work_dir / "graphite_triaxial_stress_control_case_summary.csv"
    if case_summaries:
        write_case_summary_csv(case_summary_csv, case_summaries)

    commands_csv = work_dir / "graphite_triaxial_sweep_cases.csv"
    with commands_csv.open("w", newline="") as stream:
        fieldnames = ["thetaDeg", "pressureGPa", "controlPressureGPa", "stressControlMinStrain", "stressControlMaxStrain", "outputCsv", "diagnosticCsv", "returncode", "command"]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for command in commands:
            writer.writerow(command)

    summary: Dict[str, Any] = {
        "name": "graphite_triaxial_sweep_material_point",
        "material": args.material_name,
        "initialDensity": initial_density,
        "materialPressureElasticityWarning": material_pressure_warning,
        "thetaDeg": [float(theta) for theta in theta_values],
        "pressureGPa": [float(pressure) for pressure in pressure_values],
        "zeroPressureControl": args.zero_pressure_control,
        "ambientPressureGPa": float(args.ambient_pressure_gpa),
        "regularizedZeroPressureGPa": float(args.regularized_zero_pressure_gpa),
        "dt": float(args.dt),
        "nSteps": int(args.n_steps),
        "strainRatePerMicrosecond": float(args.strain_rate),
        "targetFinalEpsZz": target_final_ezz,
        "temperatureK": 300.0,
        "controlProtocol": control_protocol_name(args.transverse_control),
        "transverseControl": args.transverse_control,
        "stressControlledComponents": controlled_components,
        "stressControlMinStrain": bound_for_metadata(stress_control_min_strain),
        "stressControlMaxStrain": bound_for_metadata(stress_control_max_strain),
        "driverStressFailurePolicy": args.driver_stress_failure_policy,
        "driverStressControlAlgorithm": args.driver_stress_control_algorithm,
        "driverStressControlDiagnostics": not args.no_driver_stress_control_diagnostics,
        "driverStressControlDiagnosticsLevel": args.driver_stress_control_diagnostics_level,
        "stressDiagnosticComponents": list(args.stress_diagnostic_components),
        "driver": args.driver,
        "preparedCases": prepared_cases,
        "caseSummaries": case_summaries,
        "caseSummaryCsv": str(case_summary_csv) if case_summaries else None,
        "commandsCsv": str(commands_csv),
        "diagnosticsCsv": str(diagnostics_csv) if diagnostic_rows else None,
        "diagnosticsSummaryCsv": str(diagnostics_summary_csv) if diagnostic_rows else None,
        "processedCsv": str(processed_csv) if processed_rows else None,
        "plots": plots,
        "passed": (not failures) and (not args.prepare_only),
        "preparedOnly": bool(args.prepare_only),
        "failures": failures,
        "warnings": warnings,
    }
    if args.prepare_only:
        summary["passed"] = True

    summary_path = work_dir / "graphite_triaxial_sweep_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", default=default_compiled_driver_executable())
    parser.add_argument("--work-dir", default="graphite_triaxial_sweep_material_point")
    parser.add_argument("--material-name", default=DEFAULT_MATERIAL_NAME)
    parser.add_argument("--theta-deg", type=parse_float_list, default=list(DEFAULT_THETA_DEG),
                        help="comma-separated material orientation angles in degrees; default: 0,30,45,60,90")
    parser.add_argument("--pressure-gpa", type=parse_float_list, default=list(DEFAULT_PRESSURE_GPA),
                        help="comma-separated confining pressures in GPa; default: 0,1,3,5,10,20,30")
    parser.add_argument("--zero-pressure-control", default="auto", choices=("auto", "exact", "ambient", "regularized"),
                        help="stress target used for the case labelled P0=0: auto uses ambient pressure for plane-strain and --regularized-zero-pressure-gpa for full-stress; exact uses 0 GPa; ambient uses 1 atm by default")
    parser.add_argument("--ambient-pressure-gpa", type=float, default=DEFAULT_AMBIENT_PRESSURE_GPA,
                        help="ambient pressure in GPa used when --zero-pressure-control ambient; default is 1 atm")
    parser.add_argument("--regularized-zero-pressure-gpa", type=float, default=1.0e-3,
                        help="small positive pressure target in GPa used when --zero-pressure-control regularized")
    parser.add_argument("--dt", type=float, default=5.0e-7, help="time step in microseconds")
    parser.add_argument("--n-steps", type=int, default=1000)
    parser.add_argument("--strain-rate", type=float, default=-1.0e3, help="axial true strain rate in 1/microsecond")
    parser.add_argument("--transverse-control", default="plane-strain", type=normalize_transverse_control,
                        help="transverse loading protocol: plane-strain keeps eps_xx=0 and stress-controls yy; full-stress stress-controls both xx and yy")
    parser.add_argument("--stress-control-min-strain", default=None,
                        help="stress-controlled total strain lower bound; use component syntax such as xx=0,yy=0 or six Voigt values. Default is xx=0,yy=0 for full-stress triaxial compression and disabled otherwise")
    parser.add_argument("--stress-control-max-strain", default=None,
                        help="stress-controlled total strain upper bound; use component syntax such as xx=0.5,yy=0.5 or six Voigt values")
    parser.add_argument("--initial-density", type=float, default=None,
                        help="initial density; default uses pfw_materials.py defaultDensity for --material-name")
    parser.add_argument("--length-scale", type=float, default=1.0e-4)
    parser.add_argument("--strength-scale", type=float, default=1.0)
    parser.add_argument("--driver-stress-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--fd-epsilon", type=float, default=1.0e-6)
    parser.add_argument("--max-newton-iterations", type=int, default=40)
    parser.add_argument("--max-line-search-iterations", type=int, default=16)
    parser.add_argument("--max-stress-bracket-iterations", type=int, default=32)
    parser.add_argument("--max-stress-bisection-iterations", type=int, default=64)
    parser.add_argument("--stress-bracket-initial-scale", type=float, default=0.0,
                        help="initial scalar stress-control bracket half-width; 0 uses the current prescribed strain scale")
    parser.add_argument("--stress-bracket-max-strain", type=float, default=5.0e-2,
                        help="maximum per-step scalar stress-control strain search distance")
    parser.add_argument("--stress-bracket-growth", type=float, default=2.0)
    parser.add_argument("--driver-stress-control-algorithm", default="hybrid",
                        choices=("newton", "regularizedNewton", "servo", "hybrid"),
                        help="stress-control update used by the compiled driver; hybrid tries Newton/LM, scalar bracketing, then a bounded servo")
    parser.add_argument("--stress-control-regularization", type=float, default=1.0e-12,
                        help="relative Levenberg-Marquardt regularization for stress-control Jacobians")
    parser.add_argument("--stress-control-max-strain-correction", type=float, default=2.0e-2,
                        help="maximum strain correction per stress-control unknown in one load step")
    parser.add_argument("--stress-control-servo-compliance", type=float, default=1.0e-2,
                        help="fallback servo compliance in strain/GPa")
    parser.add_argument("--stress-control-servo-relaxation", type=float, default=0.5,
                        help="fallback servo relaxation factor")
    parser.add_argument("--stress-control-servo-derivative-floor", type=float, default=1.0e-8,
                        help="minimum usable finite-difference stress derivative magnitude for servo fallback")
    parser.add_argument("--stress-control-servo-iterations", type=int, default=24,
                        help="maximum fallback servo attempts per load step")
    parser.add_argument("--stress-control-pattern-iterations", type=int, default=24,
                        help="derivative-free pattern-search fallback iterations per load step; useful for multi-component stress control")
    parser.add_argument("--stress-control-pattern-initial-step", type=float, default=0.0,
                        help="initial pattern-search strain step; 0 uses the current prescribed strain scale")
    parser.add_argument("--stress-control-pattern-min-step", type=float, default=1.0e-12,
                        help="minimum pattern-search strain step before giving up")
    parser.add_argument("--stress-control-pattern-shrink", type=float, default=0.5,
                        help="pattern-search shrink factor after no improvement")
    parser.add_argument("--stress-control-pattern-growth", type=float, default=1.25,
                        help="pattern-search growth factor after successful improvement")
    parser.add_argument("--driver-stress-failure-policy", default="continue", choices=("error", "stop", "continue"),
                        help="compiled-driver response to stress-control failure; default keeps writing a full curve with non-converged rows flagged")
    parser.add_argument("--no-driver-stress-control-diagnostics", action="store_true",
                        help="do not ask the compiled driver to write per-trial stress-control diagnostics")
    parser.add_argument("--driver-stress-control-diagnostics-level", default="iteration", choices=("step", "iteration", "full"),
                        help="verbosity of compiled-driver stress-control diagnostics; full includes finite-difference perturbations")
    parser.add_argument("--stress-control-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--plane-strain-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--final-strain-tolerance", type=float, default=5.0e-10)
    parser.add_argument("--prepare-only", action="store_true", help="write sidecar files but do not run the compiled driver")
    parser.add_argument("--postprocess-only", action="store_true", help="read existing case CSV files and regenerate processed CSV/plots")
    parser.add_argument("--reuse-existing", action="store_true", help="skip a driver run when the expected case output CSV already exists")
    parser.add_argument("--keep-going", action="store_true", default=True, help="continue the sweep after a case failure; enabled by default")
    parser.add_argument("--fail-fast", action="store_true", help="raise immediately on the first failed or missing case")
    parser.add_argument("--strict", action="store_true", help="return a nonzero exit code if any case is partial or fails a sanity check")
    parser.add_argument("--no-check", action="store_true", help="deprecated alias for the default non-strict behavior")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--linear-y", action="store_true", help="use a linear y-axis instead of the Figure-5-style logarithmic axis")
    parser.add_argument("--y-floor", type=float, default=1.0e-8, help="smallest positive plotted y value for log-scale plots")
    parser.add_argument("--stress-diagnostic-components", default="controlled",
                        help="comma-separated stress components for five-panel stress-vs-time diagnostics; use 'all' for all six or 'controlled' for the current stress-controlled set")
    args = parser.parse_args()

    summary = run_sweep(args)
    work_dir = Path(args.work_dir)
    print("summary JSON:", work_dir / "graphite_triaxial_sweep_summary.json")
    print("case command CSV:", summary["commandsCsv"])
    if summary.get("caseSummaryCsv"):
        print("case summary CSV:", summary["caseSummaryCsv"])
    if summary.get("diagnosticsCsv"):
        print("stress-control diagnostics CSV:", summary["diagnosticsCsv"])
    if summary.get("processedCsv"):
        print("processed CSV:", summary["processedCsv"])
    for label, path in summary.get("plots", {}).items():
        print(f"plot {label}:", path)
    print("passed:", summary["passed"])
    if summary.get("warnings"):
        for warning in summary["warnings"]:
            print("warning:", warning)
    if summary.get("failures"):
        for failure in summary["failures"]:
            print("failure:", failure)
    return 0 if summary["passed"] or args.no_check or not args.strict else 1


if __name__ == "__main__":
    raise SystemExit(main())
