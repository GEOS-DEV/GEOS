#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_simple_post import compact_float, latex_escape, read_csv_numeric, write_json, write_rows, plot_metric

DEFAULT_CASE_ID = "weakTraceRotatedBar"
VARIANT_LABELS = {
    "singleField": "single-field baseline",
    "falseElasticCZ": "false-elastic CZ reference",
    "traceContactGroups": "trace projection, contact groups",
    "dfgSurfaceFlags": "DFG surface flags placeholder",
    "stressResidualPartition": "stress-residual partition placeholder",
}
LAYER_ORDER = ["compliant_near", "stiff_near", "compliant_far", "stiff_far"]
TRACERS_PER_LAYER = 3
N = (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), 0.0)
T = (-1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), 0.0)
COHESIVE_FLAG = 3
DAMAGED_COHESIVE_FLAG = 4
WEAK_DISCONTINUITY_FLAG = 5
INTERFACE_SURFACE_FLAGS = {COHESIVE_FLAG, DAMAGED_COHESIVE_FLAG, WEAK_DISCONTINUITY_FLAG}
SURFACE_FLAG_LABELS = {
    COHESIVE_FLAG: "cohesive",
    DAMAGED_COHESIVE_FLAG: "damaged_cohesive",
    WEAK_DISCONTINUITY_FLAG: "weak_discontinuity",
}
SERIES_S0 = -0.82
SERIES_S1 = 0.82
SERIES_INTERFACE = 0.006875
DEFAULT_FINAL_STRAIN = 0.006
DEFAULT_COMPLIANT = {"bulk": 2.0 / 3.0, "shear": 0.4, "yield": 1.0}
DEFAULT_STIFF = {"bulk": 20.0 / 3.0, "shear": 4.0, "yield": 0.025}
NODE_RESTART_ROOT = "Problem/domain/MeshBodies/backgroundGrid/meshLevels/Level0/nodeManager"
INTERFACE_MECHANISM_LABELS = {
    0: "none",
    1: "cohesive_zone",
    2: "ordinary_contact",
    3: "weak_trace_anchor",
    4: "weak_trace_support",
    5: "contact_suppressed_trace_inactive",
}
NODAL_INTERFACE_FIELDS = [
    "gridInterfaceNormalForce",
    "gridInterfaceTangentialForce",
    "gridInterfaceNormalVelocityJump",
    "gridInterfaceTangentialVelocityJump",
    "gridInterfaceNormalDisplacementJump",
    "gridInterfaceTangentialDisplacementJump",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Post-process weakTraceRotatedBar verification")
    p.add_argument("--source-dir", default=str(SOURCE_DIR))
    p.add_argument("--output-dir", default=str(SOURCE_DIR / "output" / DEFAULT_CASE_ID))
    p.add_argument("--case-id", default=DEFAULT_CASE_ID)
    p.add_argument("--output-prefix", default="weakTraceRotatedBar")
    p.add_argument("--run-dir", default="")
    p.add_argument("--python", dest="python_cmd", default=sys.executable)
    p.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    p.add_argument("--visit-timeout", type=float, default=float(os.environ.get("MPM_VV_VISIT_TIMEOUT", "90")))
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def find_manifest(source_dir: Path, output_dir: Path, case_id: str) -> dict:
    for path in [output_dir / "weakTraceRotatedBar_jobs.json", output_dir / f"{case_id}_jobs.json"]:
        if path.is_file():
            return json.loads(path.read_text())
    subcases = []
    for name, label in VARIANT_LABELS.items():
        jobs_path = source_dir / "output" / f"{case_id}__{name}" / f"weakTraceRotatedBar_{name}_jobs.json"
        jobs = {}
        if jobs_path.is_file():
            try:
                jobs = json.loads(jobs_path.read_text())
            except Exception:
                jobs = {}
        subcases.append({"name": name, "label": label, "case_id": f"{case_id}__{name}", "case_name": f"weakTraceRotatedBar_{name}", "jobs": jobs, "placeholder": name not in {"singleField", "falseElasticCZ", "traceContactGroups"}})
    return {"case_id": case_id, "subcases": subcases}


def run_dir_from_subcase(source_dir: Path, subcase: dict) -> Path | None:
    run_dir = str(subcase.get("jobs", {}).get("run_dir", "")).strip()
    if run_dir and Path(run_dir).is_dir():
        return Path(run_dir)
    for p in [source_dir / "output" / subcase.get("case_id", ""), source_dir / "output" / subcase.get("case_id", "") / subcase.get("case_name", "")]:
        if p.is_dir():
            return p
    return None


def tracer_files(run_dir: Path) -> list[Path]:
    return sorted(set(run_dir.glob("weakTraceRotatedBarTracer_*.csv")) | set(run_dir.glob("tracerPoint_*.csv")))


def parse_tracer(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            out = {}
            for k, v in row.items():
                if v in (None, ""):
                    continue
                try:
                    out[k] = float(v)
                except Exception:
                    out[k] = v
            rows.append(out)
    return rows


def stress_from_row(row: dict) -> tuple[float, float, float, float, float, float] | None:
    keys = ["stressXX", "stressYY", "stressZZ", "stressYZ", "stressXZ", "stressXY"]
    if not all(k in row for k in keys):
        return None
    return tuple(float(row[k]) for k in keys)


def von_mises(sig) -> float:
    sxx, syy, szz, syz, sxz, sxy = sig
    return math.sqrt(max(0.0, 0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2) + 3.0 * (sxy*sxy + sxz*sxz + syz*syz)))


def normal_stress(sig) -> float:
    sxx, syy, szz, syz, sxz, sxy = sig
    nx, ny, nz = N
    return nx*nx*sxx + ny*ny*syy + nz*nz*szz + 2.0*(nx*ny*sxy + nx*nz*sxz + ny*nz*syz)


def tracer_layer(index: int) -> str:
    return LAYER_ORDER[min(len(LAYER_ORDER) - 1, index // TRACERS_PER_LAYER)]


def compact_int(value: object) -> str:
    try:
        return str(int(value))
    except Exception:
        return "--"


def dot3(a, b) -> float:
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]) + float(a[2]) * float(b[2])


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text or "")]


def elastic_from_bulk_shear(bulk: float, shear: float) -> tuple[float, float]:
    young = 9.0 * bulk * shear / (3.0 * bulk + shear)
    poisson = (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear))
    return young, poisson


def plane_strain_axial_modulus(bulk: float, shear: float) -> float:
    young, poisson = elastic_from_bulk_shear(bulk, shear)
    return young * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson))


def run_xml_path(run_dir: Path | None) -> Path | None:
    if run_dir is None:
        return None
    for candidate in [run_dir / "mpm_weakTraceRotatedBar.xml", *sorted(run_dir.glob("mpm_*.xml"))]:
        if candidate.is_file():
            return candidate
    return None


def analytical_solution(run_dir: Path | None = None) -> dict:
    final_strain = DEFAULT_FINAL_STRAIN
    strain_table_times = [0.0, 20.0]
    strain_table_aa = [0.0, DEFAULT_FINAL_STRAIN]
    strain_interp = "Cosine"
    compliant = dict(DEFAULT_COMPLIANT)
    stiff = dict(DEFAULT_STIFF)
    xml_path = run_xml_path(run_dir)
    if xml_path is not None:
        try:
            root = ET.parse(xml_path).getroot()
            solver = root.find(".//SolidMechanics_MPM")
            if solver is not None:
                vals = parse_floats(solver.attrib.get("fTable", ""))
                if len(vals) >= 8:
                    row_width = 4
                    rows = [vals[i:i + row_width] for i in range(0, len(vals), row_width) if len(vals[i:i + row_width]) == row_width]
                    strain_table_times = [r[0] for r in rows]
                    strain_table_aa = [0.5 * ((r[1] - 1.0) + (r[2] - 1.0)) for r in rows]
                    final_row = vals[-row_width:]
                    final_strain = 0.5 * ((final_row[1] - 1.0) + (final_row[2] - 1.0))
                strain_interp = solver.attrib.get("fTableInterpType", strain_interp)
            for mat in root.findall(".//VonMisesJ"):
                target = None
                if mat.attrib.get("name") == "weakTraceCompliantVonMises":
                    target = compliant
                elif mat.attrib.get("name") == "weakTraceStiffWeakVonMises":
                    target = stiff
                if target is not None:
                    target["bulk"] = float(mat.attrib.get("defaultBulkModulus", target["bulk"]))
                    target["shear"] = float(mat.attrib.get("defaultShearModulus", target["shear"]))
                    target["yield"] = float(mat.attrib.get("defaultYieldStrength", target["yield"]))
        except Exception:
            pass

    for mat in [compliant, stiff]:
        young, poisson = elastic_from_bulk_shear(mat["bulk"], mat["shear"])
        mat["young"] = young
        mat["poisson"] = poisson
        mat["planeStrainAxialModulus"] = plane_strain_axial_modulus(mat["bulk"], mat["shear"])

    compliant_length = SERIES_INTERFACE - SERIES_S0
    stiff_length = SERIES_S1 - SERIES_INTERFACE
    total_length = SERIES_S1 - SERIES_S0
    compliant_fraction = compliant_length / total_length
    stiff_fraction = stiff_length / total_length
    sigma = final_strain / (
        compliant_fraction / compliant["planeStrainAxialModulus"] +
        stiff_fraction / stiff["planeStrainAxialModulus"]
    )
    compliant_strain = sigma / compliant["planeStrainAxialModulus"]
    stiff_strain = sigma / stiff["planeStrainAxialModulus"]
    return {
        "assumption": "bonded 1D series bar using plane-strain axial modulus E(1-nu)/((1+nu)(1-2nu))",
        "finalStrainAA": final_strain,
        "strainTableTimes": strain_table_times,
        "strainTableAA": strain_table_aa,
        "strainInterpolation": strain_interp,
        "s0": SERIES_S0,
        "s1": SERIES_S1,
        "interfaceS": SERIES_INTERFACE,
        "compliantLengthFraction": compliant_fraction,
        "stiffLengthFraction": stiff_fraction,
        "compliantYoungModulus": compliant["young"],
        "stiffYoungModulus": stiff["young"],
        "compliantPoissonRatio": compliant["poisson"],
        "stiffPoissonRatio": stiff["poisson"],
        "compliantPlaneStrainAxialModulus": compliant["planeStrainAxialModulus"],
        "stiffPlaneStrainAxialModulus": stiff["planeStrainAxialModulus"],
        "compliantYieldStrength": compliant["yield"],
        "stiffYieldStrength": stiff["yield"],
        "expectedSigmaAA": sigma,
        "expectedSigmaNN": sigma,
        "expectedCompliantStrainAA": compliant_strain,
        "expectedStiffStrainAA": stiff_strain,
        "expectedCompliantYieldRatio": abs(sigma) / compliant["yield"] if compliant["yield"] else math.nan,
        "expectedStiffYieldRatio": abs(sigma) / stiff["yield"] if stiff["yield"] else math.nan,
    }


def analytical_strain_at_time(time: float, analytic: dict) -> float:
    times = [float(v) for v in analytic.get("strainTableTimes", [0.0, 20.0])]
    strains = [float(v) for v in analytic.get("strainTableAA", [0.0, analytic.get("finalStrainAA", DEFAULT_FINAL_STRAIN)])]
    if not times or not strains or len(times) != len(strains):
        return float(analytic.get("finalStrainAA", DEFAULT_FINAL_STRAIN))
    if time <= times[0]:
        return strains[0]
    if time >= times[-1]:
        return strains[-1]
    for i in range(len(times) - 1):
        t0, t1 = times[i], times[i + 1]
        if t0 <= time <= t1:
            xi = (time - t0) / (t1 - t0) if t1 != t0 else 1.0
            if str(analytic.get("strainInterpolation", "")).lower() == "cosine":
                xi = 0.5 * (1.0 - math.cos(math.pi * xi))
            return strains[i] + xi * (strains[i + 1] - strains[i])
    return strains[-1]


def analytical_sigma_at_time(time: float, analytic: dict) -> float:
    final_strain = float(analytic.get("finalStrainAA", DEFAULT_FINAL_STRAIN))
    if final_strain == 0.0:
        return math.nan
    return float(analytic["expectedSigmaNN"]) * analytical_strain_at_time(time, analytic) / final_strain


def analytical_time_for_visit_state(state: str, analytic: dict) -> float | None:
    times = [float(v) for v in analytic.get("strainTableTimes", [0.0, 20.0])]
    if not times:
        return None
    try:
        if state == "initial":
            return times[0]
        if state == "final":
            return times[-1]
        if state == "middle":
            return 0.5 * (times[0] + times[-1])
        if state.startswith("fraction:"):
            fraction = min(1.0, max(0.0, float(state.split(":", 1)[1])))
            return times[0] + fraction * (times[-1] - times[0])
    except Exception:
        return None
    return None


def visit_color_range(quantity: str, state: str, analytic: dict) -> tuple[float | None, float | None]:
    if quantity != "sigmaNN":
        return None, None
    time = analytical_time_for_visit_state(state, analytic)
    if time is None:
        return None, None
    expected = analytical_sigma_at_time(time, analytic)
    if not math.isfinite(expected) or expected == 0.0:
        return None, None
    return min(0.0, expected), max(0.0, expected)


def layer_material(layer: str) -> str:
    return "stiff" if layer.startswith("stiff") else "compliant"


def layer_yield_strength(layer: str, analytic: dict) -> float:
    return float(analytic["stiffYieldStrength"] if layer_material(layer) == "stiff" else analytic["compliantYieldStrength"])


def add_analytical_layer_fields(row: dict, analytic: dict) -> None:
    expected = analytical_sigma_at_time(float(row.get("time", 0.0)), analytic)
    row["expectedSigmaNN"] = expected
    row["sigmaNNError"] = float(row["meanSigmaNN"]) - expected
    row["sigmaNNRelativeError"] = row["sigmaNNError"] / expected if expected != 0.0 else math.nan
    row["sigmaNNMeasuredOverExpected"] = float(row["meanSigmaNN"]) / expected if expected != 0.0 else math.nan
    yield_strength = layer_yield_strength(str(row.get("layer", "")), analytic)
    row["yieldStrength"] = yield_strength
    row["vonMisesYieldRatio"] = float(row["meanVonMisesStress"]) / yield_strength if yield_strength else math.nan


def read_restart_array(handle, base: str, name: str):
    group = handle[f"{base}/{name}"]
    dims = tuple(int(v) for v in group["__dimensions__"][...])
    return group["__values__"][...].reshape(dims)


def latest_restart_dir(run_dir: Path) -> Path | None:
    restart_dirs = [p for p in sorted(run_dir.glob("mpm_weakTraceRotatedBar_restart_*")) if p.is_dir()]
    return restart_dirs[-1] if restart_dirs else None


def capitalized_metric_name(key: str) -> str:
    return key[0].upper() + key[1:] if key else key


def read_restart_node_array(restart_dir: Path, name: str):
    try:
        import numpy as np
        import h5py
    except Exception as exc:
        raise RuntimeError(f"h5py/numpy unavailable: {exc}") from exc

    arrays = []
    for rank_path in sorted(restart_dir.glob("rank_*.hdf5")):
        with h5py.File(rank_path, "r") as handle:
            path = f"{NODE_RESTART_ROOT}/{name}"
            if path not in handle:
                continue
            arrays.append(read_restart_array(handle, NODE_RESTART_ROOT, name))
    if not arrays:
        return None
    values = np.concatenate(arrays, axis=0)
    while len(values.shape) > 1 and values.shape[-1] == 1:
        values = values[..., 0]
    return values


def percentile_abs(values, percentile: float) -> float:
    vals = sorted(abs(float(v)) for v in values if math.isfinite(float(v)))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * percentile / 100.0
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    return vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)


def process_nodal_interface_diagnostics(subcase: dict, run_dir: Path | None) -> tuple[list[dict], dict]:
    summary = {
        "nodalInterfaceStatus": "not evaluated",
    }
    if run_dir is None:
        summary["nodalInterfaceStatus"] = "missing run directory"
        return [], summary
    restart_dir = latest_restart_dir(run_dir)
    if restart_dir is None:
        summary["nodalInterfaceStatus"] = "missing restart directory"
        return [], summary

    try:
        mechanism = read_restart_node_array(restart_dir, "gridInterfaceMechanism")
    except Exception as exc:
        summary["nodalInterfaceStatus"] = str(exc)
        return [], summary
    if mechanism is None:
        summary["nodalInterfaceStatus"] = "missing gridInterfaceMechanism in final restart"
        return [], summary

    fields = {}
    for name in NODAL_INTERFACE_FIELDS:
        try:
            fields[name] = read_restart_node_array(restart_dir, name)
        except Exception:
            fields[name] = None

    rows = []
    num_fields = int(mechanism.shape[1]) if len(mechanism.shape) > 1 else 1
    aggregate_counts = {code: 0 for code in INTERFACE_MECHANISM_LABELS}
    aggregate_interface_entries = 0
    aggregate_values = {name: [] for name in NODAL_INTERFACE_FIELDS}

    for field_index in range(num_fields):
        mech = mechanism[:, field_index] if num_fields > 1 else mechanism[:]
        active_mask = [int(v) != 0 for v in mech]
        row = {
            "variant": subcase.get("name"),
            "label": subcase.get("label"),
            "fieldIndex": field_index,
            "fieldLabel": f"velocityField{field_index + 1}",
            "restartDir": str(restart_dir),
            "interfaceEntries": int(sum(active_mask)),
        }
        for code, label in INTERFACE_MECHANISM_LABELS.items():
            count = int(sum(1 for v in mech if int(v) == code))
            row[f"mechanism_{code}_{label}_count"] = count
            aggregate_counts[code] += count
        aggregate_interface_entries += row["interfaceEntries"]

        for name, values in fields.items():
            if values is None:
                row[f"{name}_status"] = "missing"
                continue
            vals = values[:, field_index] if len(values.shape) > 1 else values[:]
            active_vals = [float(v) for v, active in zip(vals, active_mask) if active]
            aggregate_values[name].extend(active_vals)
            row[f"{name}_meanAbs"] = sum(abs(v) for v in active_vals) / len(active_vals) if active_vals else math.nan
            row[f"{name}_p95Abs"] = percentile_abs(active_vals, 95.0)
            row[f"{name}_maxAbs"] = max((abs(v) for v in active_vals), default=math.nan)
        rows.append(row)

    summary.update({
        "nodalInterfaceStatus": "processed",
        "nodalInterfaceRestart": str(restart_dir),
        "nodalInterfaceFieldCount": num_fields,
        "nodalInterfaceEntries": aggregate_interface_entries,
    })
    for code, label in INTERFACE_MECHANISM_LABELS.items():
        summary[f"nodalMechanism_{code}_{label}_count"] = aggregate_counts[code]
    for name, vals in aggregate_values.items():
        summary[f"nodal_{name}_p95Abs"] = percentile_abs(vals, 95.0)
        summary[f"nodal_{name}_maxAbs"] = max((abs(v) for v in vals), default=math.nan)
    return rows, summary


def restart_surface_points(restart_dir: Path, surface_flags: set[int] | None = None) -> list[dict]:
    try:
        import h5py
    except Exception as exc:
        raise RuntimeError(f"h5py unavailable: {exc}") from exc

    accepted_flags = set(INTERFACE_SURFACE_FLAGS if surface_flags is None else surface_flags)
    points = []
    root = "Problem/domain/MeshBodies/particles/meshLevels/Level0/ParticleRegions/particleRegionsGroup"
    for rank_path in sorted(restart_dir.glob("rank_*.hdf5")):
        with h5py.File(rank_path, "r") as handle:
            if root not in handle:
                continue
            for region_name in handle[root]:
                region_path = f"{root}/{region_name}"
                region_group = handle[region_path]
                if not isinstance(region_group, h5py.Group) or "particleSubRegions" not in region_group:
                    continue
                subregions = handle[f"{region_path}/particleSubRegions"]
                for sub_name in subregions:
                    base = f"{region_path}/particleSubRegions/{sub_name}"
                    sub_group = handle[base]
                    if not isinstance(sub_group, h5py.Group) or "particleSurfaceFlag" not in sub_group:
                        continue
                    center = read_restart_array(handle, base, "particleCenter")
                    surface_position = read_restart_array(handle, base, "particleSurfacePosition")
                    surface_flag = read_restart_array(handle, base, "particleSurfaceFlag")
                    particle_group = read_restart_array(handle, base, "particleGroup")
                    material_type = read_restart_array(handle, base, "particleMaterialType")
                    for p, flag in enumerate(surface_flag):
                        flag_int = int(flag)
                        if flag_int not in accepted_flags:
                            continue
                        point = center[p] + surface_position[p]
                        material = int(material_type[p])
                        group = int(particle_group[p])
                        side = material if material in (0, 1) else group
                        points.append({
                            "key": f"{rank_path.name}:{region_name}:{sub_name}:{p}",
                            "rank": rank_path.name,
                            "region": region_name,
                            "subregion": sub_name,
                            "particleIndex": int(p),
                            "group": group,
                            "materialType": material,
                            "side": side,
                            "surfaceFlag": flag_int,
                            "surfaceFlagLabel": SURFACE_FLAG_LABELS.get(flag_int, str(flag_int)),
                            "point": [float(point[0]), float(point[1]), float(point[2])],
                        })
    return points


def match_interface_pairs(initial_points: list[dict]) -> list[tuple[dict, dict, float]]:
    pairs = []
    for flag in sorted({int(p["surfaceFlag"]) for p in initial_points}):
        lower = [p for p in initial_points if int(p["surfaceFlag"]) == flag and int(p.get("side", -1)) == 0]
        upper = [p for p in initial_points if int(p["surfaceFlag"]) == flag and int(p.get("side", -1)) == 1]
        unused = set(range(len(upper)))
        for lo in sorted(lower, key=lambda p: (dot3(p["point"], T), float(p["point"][2]))):
            best = None
            lo_t = dot3(lo["point"], T)
            lo_z = float(lo["point"][2])
            for idx in list(unused):
                hi = upper[idx]
                dt = dot3(hi["point"], T) - lo_t
                dz = float(hi["point"][2]) - lo_z
                dist2 = dt * dt + dz * dz
                if best is None or dist2 < best[0]:
                    best = (dist2, idx)
            if best is not None:
                unused.remove(best[1])
                pairs.append((lo, upper[best[1]], math.sqrt(best[0])))
    return pairs


def solver_attributes(run_dir: Path | None) -> dict:
    xml_path = run_xml_path(run_dir)
    if xml_path is None:
        return {}
    try:
        root = ET.parse(xml_path).getroot()
        solver = root.find(".//SolidMechanics_MPM")
        if solver is None:
            return {}
        keys = [
            "explicitSurfaceNormalInfluence",
            "weakInterfaceTraceProjectionScale",
            "weakInterfaceTraceProjectionIterations",
            "weakInterfaceTraceGapStabilization",
            "weakInterfaceTraceMaxGapCorrectionVelocity",
            "weakInterfaceTraceSuppressNodalContact",
        ]
        return {f"solver_{key}": solver.attrib.get(key, "") for key in keys if key in solver.attrib}
    except Exception:
        return {}


def process_interface_kinematics(subcase: dict, run_dir: Path | None) -> tuple[list[dict], dict]:
    summary = {
        "interfaceKinematicsStatus": "not evaluated",
        "interfaceMatchedPairs": 0,
    }
    if run_dir is None:
        summary["interfaceKinematicsStatus"] = "missing run directory"
        return [], summary
    restart_dirs = [p for p in sorted(run_dir.glob("mpm_weakTraceRotatedBar_restart_*")) if p.is_dir()]
    if len(restart_dirs) < 2:
        summary["interfaceKinematicsStatus"] = "missing initial/final restart directories"
        return [], summary

    try:
        initial_points = restart_surface_points(restart_dirs[0])
        final_points = restart_surface_points(restart_dirs[-1])
    except Exception as exc:
        summary["interfaceKinematicsStatus"] = str(exc)
        return [], summary

    final_by_key = {p["key"]: p for p in final_points}
    pairs = match_interface_pairs(initial_points)
    rows = []
    for pair_index, (lower, upper, match_distance) in enumerate(pairs):
        lower_final = final_by_key.get(lower["key"])
        upper_final = final_by_key.get(upper["key"])
        if lower_final is None or upper_final is None:
            continue
        initial_jump = [float(upper["point"][i]) - float(lower["point"][i]) for i in range(3)]
        final_jump = [float(upper_final["point"][i]) - float(lower_final["point"][i]) for i in range(3)]
        initial_mid = [0.5 * (float(upper["point"][i]) + float(lower["point"][i])) for i in range(3)]
        final_mid = [0.5 * (float(upper_final["point"][i]) + float(lower_final["point"][i])) for i in range(3)]
        initial_normal_gap = dot3(initial_jump, N)
        final_normal_gap = dot3(final_jump, N)
        initial_tangential_slip = dot3(initial_jump, T)
        final_tangential_slip = dot3(final_jump, T)
        rows.append({
            "variant": subcase.get("name"),
            "label": subcase.get("label"),
            "pairIndex": pair_index,
            "lowerParticle": lower["key"],
            "upperParticle": upper["key"],
            "surfaceFlag": int(lower.get("surfaceFlag", -1)),
            "surfaceFlagLabel": lower.get("surfaceFlagLabel", ""),
            "initialMatchDistance": match_distance,
            "initialInterfaceT": dot3(initial_mid, T),
            "finalInterfaceT": dot3(final_mid, T),
            "initialInterfaceN": dot3(initial_mid, N),
            "finalInterfaceN": dot3(final_mid, N),
            "initialNormalGap": initial_normal_gap,
            "finalNormalGap": final_normal_gap,
            "normalGapChange": final_normal_gap - initial_normal_gap,
            "initialTangentialSlip": initial_tangential_slip,
            "finalTangentialSlip": final_tangential_slip,
            "tangentialSlipChange": final_tangential_slip - initial_tangential_slip,
        })

    summary.update({
        "interfaceKinematicsStatus": "processed" if rows else "no matched weak-interface particles",
        "interfaceInitialParticleCount": len(initial_points),
        "interfaceFinalParticleCount": len(final_points),
        "interfaceInitialWeakParticleCount": sum(1 for p in initial_points if int(p.get("surfaceFlag", -1)) == WEAK_DISCONTINUITY_FLAG),
        "interfaceFinalWeakParticleCount": sum(1 for p in final_points if int(p.get("surfaceFlag", -1)) == WEAK_DISCONTINUITY_FLAG),
        "interfaceInitialCohesiveParticleCount": sum(1 for p in initial_points if int(p.get("surfaceFlag", -1)) in (COHESIVE_FLAG, DAMAGED_COHESIVE_FLAG)),
        "interfaceFinalCohesiveParticleCount": sum(1 for p in final_points if int(p.get("surfaceFlag", -1)) in (COHESIVE_FLAG, DAMAGED_COHESIVE_FLAG)),
        "interfaceMatchedPairs": len(rows),
        "interfaceRestartInitial": str(restart_dirs[0]),
        "interfaceRestartFinal": str(restart_dirs[-1]),
    })
    if rows:
        for key in ["initialNormalGap", "finalNormalGap", "normalGapChange", "initialTangentialSlip", "finalTangentialSlip", "tangentialSlipChange"]:
            vals = [float(r[key]) for r in rows]
            cap = capitalized_metric_name(key)
            summary[f"interfaceMean{cap}"] = sum(vals) / len(vals)
            summary[f"interfaceMin{cap}"] = min(vals)
            summary[f"interfaceMax{cap}"] = max(vals)
            summary[f"interfaceMeanAbs{cap}"] = sum(abs(v) for v in vals) / len(vals)
            summary[f"interfaceMaxAbs{cap}"] = max(abs(v) for v in vals)
            max_row = max(rows, key=lambda r: abs(float(r[key])))
            summary[f"interfaceMaxAbs{cap}PairIndex"] = int(max_row.get("pairIndex", -1))
            summary[f"interfaceMaxAbs{cap}T"] = float(max_row.get("initialInterfaceT", math.nan))
    return rows, summary


def process_variant(source_dir: Path, subcase: dict, analytic: dict) -> tuple[list[dict], list[dict], list[dict], list[dict], dict]:
    if subcase.get("placeholder"):
        return [], [], [], [], {"variant": subcase.get("name"), "label": subcase.get("label"), "status": "placeholder"}
    run_dir = run_dir_from_subcase(source_dir, subcase)
    if run_dir is None:
        return [], [], [], [], {"variant": subcase.get("name"), "label": subcase.get("label"), "status": "missing run directory"}
    interface_rows, interface_summary = process_interface_kinematics(subcase, run_dir)
    nodal_rows, nodal_summary = process_nodal_interface_diagnostics(subcase, run_dir)
    files = tracer_files(run_dir)
    if not files:
        summary = {"variant": subcase.get("name"), "label": subcase.get("label"), "status": "missing tracer files", "run_dir": str(run_dir)}
        summary.update(interface_summary)
        summary.update(nodal_summary)
        return [], [], interface_rows, nodal_rows, summary

    rows = []
    for idx, path in enumerate(files):
        layer = tracer_layer(idx)
        for row in parse_tracer(path):
            sig = stress_from_row(row)
            if sig is None:
                continue
            t = float(row.get("t", row.get("time", 0.0)))
            rows.append({
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "tracer": path.name,
                "tracerIndex": idx,
                "layer": layer,
                "time": t,
                "sigmaNN": normal_stress(sig),
                "vonMisesStress": von_mises(sig),
                "plasticStrainMagnitude": float(row.get("plasticStrainMagnitude", 0.0) or 0.0),
                "materialType": row.get("materialType", ""),
                "density": row.get("density", ""),
            })

    by = defaultdict(list)
    for r in rows:
        key = (r["variant"], r["label"], r["layer"], r["time"])
        by[key].append(r)
    layer_rows = []
    for (variant, label, layer, t), vals in sorted(by.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][3])):
        n = len(vals)
        layer_row = {
            "variant": variant,
            "label": label,
            "layer": layer,
            "time": t,
            "numTracers": n,
            "meanSigmaNN": sum(v["sigmaNN"] for v in vals) / n,
            "meanVonMisesStress": sum(v["vonMisesStress"] for v in vals) / n,
            "meanPlasticStrainMagnitude": sum(v["plasticStrainMagnitude"] for v in vals) / n,
        }
        add_analytical_layer_fields(layer_row, analytic)
        layer_rows.append(layer_row)
    final_by_layer = {r["layer"]: r for r in layer_rows if r["time"] == max((x["time"] for x in layer_rows), default=0.0)}
    summary = {
        "variant": subcase.get("name"),
        "label": subcase.get("label"),
        "status": "processed",
        "run_dir": str(run_dir),
        "numTracerRows": len(rows),
    }
    summary.update(solver_attributes(run_dir))
    for layer in LAYER_ORDER:
        if layer in final_by_layer:
            summary[f"final_{layer}_vm"] = final_by_layer[layer]["meanVonMisesStress"]
            summary[f"final_{layer}_sigmaNN"] = final_by_layer[layer]["meanSigmaNN"]
            summary[f"final_{layer}_sigmaNN_error"] = final_by_layer[layer]["sigmaNNError"]
            summary[f"final_{layer}_sigmaNN_relative_error"] = final_by_layer[layer]["sigmaNNRelativeError"]
            summary[f"final_{layer}_sigmaNN_measured_over_expected"] = final_by_layer[layer]["sigmaNNMeasuredOverExpected"]
            summary[f"final_{layer}_vm_yield_ratio"] = final_by_layer[layer]["vonMisesYieldRatio"]
            summary[f"final_{layer}_plastic"] = final_by_layer[layer]["meanPlasticStrainMagnitude"]
    summary.update(interface_summary)
    summary.update(nodal_summary)
    return rows, layer_rows, interface_rows, nodal_rows, summary


def plot_layers(output_dir: Path, layer_rows: list[dict], column: str, filename: str, ylabel: str, expected_series: list[tuple[float, float]] | None = None):
    series = []
    for variant in sorted({r["variant"] for r in layer_rows}):
        for layer in LAYER_ORDER:
            data = sorted([r for r in layer_rows if r["variant"] == variant and r["layer"] == layer], key=lambda r: r["time"])
            if data:
                series.append((f"{variant} / {layer}", [r["time"] for r in data], [r[column] for r in data]))
    if expected_series:
        series.append(("bonded analytical", [p[0] for p in expected_series], [p[1] for p in expected_series]))
    plot_metric(output_dir, filename, filename.replace("_", " ").replace(".png", ""), "time", ylabel, series)


def plot_interface_gaps(output_dir: Path, interface_rows: list[dict]) -> list[str]:
    outputs = []
    for column, filename, ylabel in [
        ("normalGapChange", "weak_trace_rotated_bar_interface_normal_gap.png", "final-initial normal gap"),
        ("tangentialSlipChange", "weak_trace_rotated_bar_interface_tangential_slip.png", "final-initial tangential slip"),
    ]:
        series = []
        for variant in sorted({r["variant"] for r in interface_rows}):
            data = sorted([r for r in interface_rows if r["variant"] == variant], key=lambda r: float(r.get("initialInterfaceT", 0.0)))
            if data:
                series.append((variant, [float(r["initialInterfaceT"]) for r in data], [float(r[column]) for r in data]))
        if plot_metric(output_dir, filename, filename.replace("_", " ").replace(".png", ""), "interface tangential coordinate", ylabel, series):
            outputs.append(filename)
    return outputs


def render_visit(args, subcases: list[dict], output_dir: Path, analytic: dict):
    messages = []
    if args.no_visit:
        return [], ["VisIt render skipped: --no-visit was requested."]
    visit = args.visit_cmd or os.environ.get("VISIT_COMMAND", "") or shutil.which("visit")
    if not visit:
        return [], ["VisIt render skipped: no VISIT_COMMAND was supplied and visit was not found on PATH."]
    if not Path(visit).is_file():
        resolved = shutil.which(visit)
        if not resolved:
            return [], [f"VisIt render skipped: command {visit!r} was not found on PATH."]
        visit = resolved
    renderer = Path(args.source_dir) / "visitRender_weakTraceRotatedBar.py"
    grid_renderer = Path(args.source_dir) / "visitRender_weakTraceGridDiagnostics.py"
    frames = []
    render_specs = [
        ("vm", "final", "vm_final"),
        ("sigmaNN", "fraction:0.4", "sigmaNN_trace_band"),
    ]
    grid_render_specs = [
        ("mechanism", "final", "grid_mechanism_final", "all"),
        ("normalForce", "final", "grid_normal_force_final", "all"),
        ("normalVelocityJump", "final", "grid_normal_velocity_jump_final", "all"),
        ("normalDisplacementJump", "final", "grid_normal_displacement_jump_final", "all"),
        ("surfaceNormal", "final", "grid_surface_normal_field0_final", "0"),
        ("surfaceNormal", "final", "grid_surface_normal_field1_final", "1"),
        ("explicitSurfaceNormal", "initial", "grid_explicit_surface_normal_initial", "all"),
    ]
    for subcase in subcases:
        if subcase.get("placeholder"):
            continue
        run_dir = run_dir_from_subcase(Path(args.source_dir), subcase)
        if run_dir is None or not (run_dir / "siloFiles").is_dir():
            continue
        render_script = renderer if renderer.is_file() else run_dir / "visitRender_weakTraceRotatedBar.py"
        for quantity, state, suffix in render_specs:
            out = output_dir / f"weak_trace_rotated_bar_{subcase.get('name')}_{suffix}.png"
            color_min, color_max = visit_color_range(quantity, state, analytic)
            cmd = [
                visit,
                "-nowin",
                "-cli",
                "-s",
                str(render_script),
                "--run-dir",
                str(run_dir),
                "--output",
                str(out),
                "--state",
                state,
                "--quantity",
                quantity,
            ]
            if color_min is not None:
                cmd.extend(["--color-min", str(color_min)])
            if color_max is not None:
                cmd.extend(["--color-max", str(color_max)])
            try:
                if out.is_file():
                    out.unlink()
                result = subprocess.run(cmd, cwd=run_dir, timeout=args.visit_timeout, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                if out.is_file():
                    frames.append(out)
                    scale = f" color=[{color_min}, {color_max}]" if color_min is not None and color_max is not None else ""
                    messages.append(f"Rendered {out.name} from {subcase.get('name')} quantity={quantity} state={state}{scale}.")
                else:
                    messages.append(f"Missing VisIt frame {out.name}; returncode={result.returncode}.")
                    if result.stdout:
                        messages.append(result.stdout[-4000:])
            except Exception as exc:
                messages.append(f"Failed to render {out.name}: {exc}")
        if grid_renderer.is_file():
            for quantity, state, suffix, field_index in grid_render_specs:
                out = output_dir / f"weak_trace_rotated_bar_{subcase.get('name')}_{suffix}.png"
                cmd = [
                    visit,
                    "-nowin",
                    "-cli",
                    "-s",
                    str(grid_renderer),
                    "--run-dir",
                    str(run_dir),
                    "--output",
                    str(out),
                    "--state",
                    state,
                    "--quantity",
                    quantity,
                    "--field-index",
                    field_index,
                ]
                try:
                    if out.is_file():
                        out.unlink()
                    result = subprocess.run(cmd, cwd=run_dir, timeout=args.visit_timeout, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    if out.is_file():
                        frames.append(out)
                        messages.append(f"Rendered {out.name} from {subcase.get('name')} nodal quantity={quantity} state={state}.")
                    else:
                        messages.append(f"Missing VisIt grid frame {out.name}; returncode={result.returncode}.")
                        if result.stdout:
                            messages.append(result.stdout[-4000:])
                except Exception as exc:
                    messages.append(f"Failed to render {out.name}: {exc}")
        else:
            messages.append(f"VisIt grid render skipped: {grid_renderer} was not found.")
    return frames, messages


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = find_manifest(source_dir, output_dir, args.case_id)
    analytic_run_dir = None
    for subcase in manifest.get("subcases", []):
        if subcase.get("name") == "traceContactGroups":
            analytic_run_dir = run_dir_from_subcase(source_dir, subcase)
            break
    if analytic_run_dir is None and manifest.get("subcases"):
        analytic_run_dir = run_dir_from_subcase(source_dir, manifest.get("subcases", [])[0])
    analytic = analytical_solution(analytic_run_dir)

    all_rows = []
    all_layers = []
    all_interface_rows = []
    all_nodal_rows = []
    summaries = []
    for subcase in manifest.get("subcases", []):
        rows, layers, interface_rows, nodal_rows, summary = process_variant(source_dir, subcase, analytic)
        all_rows.extend(rows)
        all_layers.extend(layers)
        all_interface_rows.extend(interface_rows)
        all_nodal_rows.extend(nodal_rows)
        summaries.append(summary)

    write_rows(output_dir / "weak_trace_rotated_bar_tracer_metrics.csv", all_rows)
    write_rows(output_dir / "weak_trace_rotated_bar_layer_metrics.csv", all_layers)
    final_layers = []
    for variant in sorted({r["variant"] for r in all_layers}):
        max_t = max((r["time"] for r in all_layers if r["variant"] == variant), default=None)
        if max_t is not None:
            final_layers.extend([r for r in all_layers if r["variant"] == variant and r["time"] == max_t])
    write_rows(output_dir / "weak_trace_rotated_bar_final_layer_metrics.csv", final_layers)
    write_rows(output_dir / "weak_trace_rotated_bar_interface_kinematics.csv", all_interface_rows)
    write_rows(output_dir / "weak_trace_rotated_bar_nodal_interface_diagnostics.csv", all_nodal_rows)
    write_rows(output_dir / "weak_trace_rotated_bar_analytical_solution.csv", [analytic])
    write_json(output_dir / "weak_trace_rotated_bar_summary.json", {"analyticalSolution": analytic, "summaries": summaries})

    expected_curve = []
    for t in sorted({float(r["time"]) for r in all_layers}):
        expected_curve.append((t, analytical_sigma_at_time(t, analytic)))
    plot_layers(output_dir, all_layers, "meanSigmaNN", "weak_trace_rotated_bar_layer_normal_stress.png", "layer mean normal stress", expected_series=expected_curve)
    plot_layers(output_dir, all_layers, "meanVonMisesStress", "weak_trace_rotated_bar_layer_vm_stress.png", "layer mean von Mises stress")
    plot_layers(output_dir, all_layers, "meanPlasticStrainMagnitude", "weak_trace_rotated_bar_layer_plastic.png", "layer mean plastic strain")
    interface_gap_plots = plot_interface_gaps(output_dir, all_interface_rows)
    frames, visit_messages = render_visit(args, manifest.get("subcases", []), output_dir, analytic)
    if visit_messages:
        (output_dir / "weak_trace_rotated_bar_visit_render_log.txt").write_text("\n".join(visit_messages) + "\n")

    tex = []
    tex.append(r"\paragraph{Local stress diagnostics.} The rotated-bar verification compares local tracer stresses on both sides of the inclined interface.  The primary quantities are the material-normal stress $\sigma_{nn}=\mathbf{n}\cdot\boldsymbol{\sigma}\mathbf{n}$ and von Mises stress reconstructed from tracer stress components.  The stiff phase is deliberately weak so that a mixed-cell stress spike can cause spurious plasticity even though the correct series response is elastic.")
    tex.append(r"\paragraph{Analytical bonded-series solution.} The analytical reference treats the two material segments as a bonded series bar along $\mathbf{n}=(1,1,0)/\sqrt{2}$, using the plane-strain axial modulus $E(1-\nu)/((1+\nu)(1-2\nu))$ and the same F-table ramp as the input deck.  This is a local stress-transfer target for the weak trace, not a global reaction-force check.")
    tex.append(r"{\scriptsize\begin{tabular}{lrrrrrr}\toprule $\epsilon_{nn}$ & $L_c/L$ & $L_s/L$ & expected $\sigma_{nn}$ & $\epsilon_c$ & $\epsilon_s$ & stiff expected $\sigma/y$ \\\midrule")
    tex.append(
        compact_float(analytic["finalStrainAA"]) + " & " +
        compact_float(analytic["compliantLengthFraction"]) + " & " +
        compact_float(analytic["stiffLengthFraction"]) + " & " +
        compact_float(analytic["expectedSigmaNN"]) + " & " +
        compact_float(analytic["expectedCompliantStrainAA"]) + " & " +
        compact_float(analytic["expectedStiffStrainAA"]) + " & " +
        compact_float(analytic["expectedStiffYieldRatio"]) + r" \\"
    )
    tex.append(r"\bottomrule\end{tabular}}")
    tex.append(r"{\scriptsize\begin{tabular}{lrrrr}\toprule Variant & stiff near VM & stiff far VM & stiff near plastic & status \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("final_stiff_near_vm")) + " & " + compact_float(s.get("final_stiff_far_vm")) + " & " + compact_float(s.get("final_stiff_near_plastic")) + " & " + latex_escape(s.get("status", "")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    if frames:
        tex.append(r"\paragraph{VisIt field diagnostics.} The report includes particle pseudocolor frames for final von Mises stress and an intermediate-time $\sigma_{nn}=0.5\sigma_{xx}+0.5\sigma_{yy}+\sigma_{xy}$ trace-band view, plus nodal grid frames for interface mechanism, normal force, normal velocity jump, and normal displacement jump.  Grid surface normals are rendered as vector glyphs for each velocity field, with the explicit CZ partitioning normal rendered when that grid field is present.  The particle $\sigma_{nn}$ frame uses the same color limits for both material particle regions, capped by the bonded analytical stress at that render time, so local overshoots saturate instead of rescaling the two sides independently.")
    elif visit_messages:
        tex.append(r"\paragraph{VisIt field diagnostics.} VisIt field frames were requested but not generated in this post-processing environment.  See \texttt{weak\_trace\_rotated\_bar\_visit\_render\_log.txt} in the case output directory for the command or expression failure.")
    if final_layers:
        tex.append(r"{\scriptsize\begin{tabular}{llrrrr}\toprule Variant & layer & measured $\sigma_{nn}$ & expected $\sigma_{nn}$ & measured/expected & VM/yield \\\midrule")
        for r in sorted(final_layers, key=lambda x: (str(x.get("variant", "")), LAYER_ORDER.index(x.get("layer", LAYER_ORDER[-1])) if x.get("layer") in LAYER_ORDER else 99)):
            tex.append(
                latex_escape(r.get("label", r.get("variant", ""))) + " & " +
                latex_escape(r.get("layer", "")) + " & " +
                compact_float(r.get("meanSigmaNN")) + " & " +
                compact_float(r.get("expectedSigmaNN")) + " & " +
                compact_float(r.get("sigmaNNMeasuredOverExpected")) + " & " +
                compact_float(r.get("vonMisesYieldRatio")) + r" \\"
            )
        tex.append(r"\bottomrule\end{tabular}}")
    tex.append(r"{\scriptsize\begin{tabular}{lrrrrrr}\toprule Variant & pairs & mean $\Delta g_n$ & mean $|\Delta g_n|$ & max $|\Delta g_n|$ & $t$ at max & mean $|\Delta s_t|$ \\\midrule")
    for s in summaries:
        tex.append(
            latex_escape(s.get("label", s.get("variant", ""))) + " & " +
            compact_int(s.get("interfaceMatchedPairs")) + " & " +
            compact_float(s.get("interfaceMeanNormalGapChange")) + " & " +
            compact_float(s.get("interfaceMeanAbsNormalGapChange")) + " & " +
            compact_float(s.get("interfaceMaxAbsNormalGapChange")) + " & " +
            compact_float(s.get("interfaceMaxAbsNormalGapChangeT")) + " & " +
            compact_float(s.get("interfaceMeanAbsTangentialSlipChange")) + r" \\"
        )
    tex.append(r"\bottomrule\end{tabular}}")
    tex.append(r"\paragraph{Nodal interface diagnostics.} The common grid diagnostics classify each velocity-field node entry by the dominant interface mechanism: 1 cohesive zone, 2 ordinary contact, 3 weak-trace anchor, 4 weak-trace support, and 5 contact suppressed without an active trace.  The force and jump columns below are absolute values over entries with a nonzero mechanism in the final restart.")
    tex.append(r"{\scriptsize\begin{tabular}{lrrrrrrrr}\toprule Variant & CZ & contact & weak anchor & weak support & suppressed & p95 $|f_n|$ & max $|f_n|$ & p95 $|[\![v_n]\!]|$ \\\midrule")
    for s in summaries:
        tex.append(
            latex_escape(s.get("label", s.get("variant", ""))) + " & " +
            compact_int(s.get("nodalMechanism_1_cohesive_zone_count")) + " & " +
            compact_int(s.get("nodalMechanism_2_ordinary_contact_count")) + " & " +
            compact_int(s.get("nodalMechanism_3_weak_trace_anchor_count")) + " & " +
            compact_int(s.get("nodalMechanism_4_weak_trace_support_count")) + " & " +
            compact_int(s.get("nodalMechanism_5_contact_suppressed_trace_inactive_count")) + " & " +
            compact_float(s.get("nodal_gridInterfaceNormalForce_p95Abs")) + " & " +
            compact_float(s.get("nodal_gridInterfaceNormalForce_maxAbs")) + " & " +
            compact_float(s.get("nodal_gridInterfaceNormalVelocityJump_p95Abs")) + r" \\"
        )
    tex.append(r"\bottomrule\end{tabular}}")
    trace_summary = next((s for s in summaries if s.get("variant") == "traceContactGroups"), None)
    if trace_summary and trace_summary.get("interfaceKinematicsStatus") == "processed":
        measured = trace_summary.get("final_stiff_near_sigmaNN", trace_summary.get("final_compliant_near_sigmaNN"))
        measured_ratio = None
        try:
            measured_ratio = float(measured) / float(analytic["expectedSigmaNN"])
        except Exception:
            measured_ratio = None
        gap_change = trace_summary.get("interfaceMeanAbsNormalGapChange")
        slip_change = trace_summary.get("interfaceMeanAbsTangentialSlipChange")
        normal_influence = trace_summary.get("solver_explicitSurfaceNormalInfluence", "--")
        max_gap_correction = trace_summary.get("solver_weakInterfaceTraceMaxGapCorrectionVelocity", "--")
        if gap_change is not None:
            if abs(float(gap_change)) > 1.0e-4:
                cz_summary = next((s for s in summaries if s.get("variant") == "falseElasticCZ"), None)
                cz_gap_text = ""
                if cz_summary and cz_summary.get("interfaceKinematicsStatus") == "processed":
                    cz_gap_text = (
                        r"  The false-elastic CZ reference gives mean $|\Delta g_n|="
                        + compact_float(cz_summary.get("interfaceMeanAbsNormalGapChange"))
                        + r"$ and max $|\Delta g_n|="
                        + compact_float(cz_summary.get("interfaceMaxAbsNormalGapChange"))
                        + r"$ in the same restart-paired diagnostic."
                    )
                tex.append(
                    r"\paragraph{Current trace-run interpretation.} The restart-based weak-interface surface points separate in the normal direction: mean $|\Delta g_n|="
                    + compact_float(gap_change)
                    + r"$ and max $|\Delta g_n|="
                    + compact_float(trace_summary.get("interfaceMaxAbsNormalGapChange"))
                    + r"$.  The tangential mismatch is much smaller, mean $|\Delta s_t|="
                    + compact_float(slip_change)
                    + r"$."
                    + cz_gap_text
                    + r"  The trace run used explicit surface-normal influence "
                    + latex_escape(normal_influence)
                    + r" and gap stabilization "
                    + latex_escape(trace_summary.get("solver_weakInterfaceTraceGapStabilization", "--"))
                    + r" with max gap-correction velocity "
                    + latex_escape(max_gap_correction)
                    + r"; with zero gap stabilization, the projection damps velocity mismatch but does not drive accumulated surface-position opening back to zero.  The final near-interface stiff $\sigma_{nn}$ is "
                    + (compact_float(measured_ratio) if measured_ratio is not None else "--")
                    + r" of the bonded-series target.  This points to weak-trace normal separation/opening in the present run, not merely a residual stress artifact at an otherwise tied interface."
                )
            else:
                tex.append(
                    r"\paragraph{Current trace-run interpretation.} The restart-based weak-interface surface points remain nearly coincident.  The trace run used explicit surface-normal influence "
                    + latex_escape(normal_influence)
                    + r" and max gap-correction velocity "
                    + latex_escape(max_gap_correction)
                    + r", so remaining error is more consistent with a stress reconstruction or mixed-cell artifact than with macroscopic trace separation."
                )
    for name in [
        "weak_trace_rotated_bar_layer_normal_stress.png",
        "weak_trace_rotated_bar_layer_vm_stress.png",
        "weak_trace_rotated_bar_layer_plastic.png",
        *interface_gap_plots,
    ]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.49\linewidth]{\CaseOutputDir/" + name + "}")
    for frame in frames:
        try:
            rel = frame.relative_to(output_dir)
        except ValueError:
            rel = frame.name
        tex.append(r"\includegraphics[width=0.49\linewidth]{\CaseOutputDir/" + str(rel).replace(os.sep, "/") + "}")
    (output_dir / "weak_trace_rotated_bar_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
