#!/usr/bin/env python3
"""Evaluate three-field contact runs and build an expected-vs-actual report.

The script has two compatible modes:

* The standard MPM verification harness calls it with ``--run-dir`` to reduce
  one case immediately after GEOS finishes.
* Running it without ``--run-dir`` discovers all eleven case manifests and
  compiles a suite-level PDF/Markdown/CSV/JSON report.

If a junction history does not already exist, the script invokes
``visitExtract_threeFieldContact.py`` with VisIt to sample the grid node nearest
the origin from the Silo time series.
"""

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
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Vec3 = Tuple[float, float, float]
BODY_NAMES = ("A", "B", "C")
COMPONENTS = ("x", "y", "z")


@dataclass(frozen=True)
class CaseSpec:
    key: str
    title: str
    edge_case: str
    family: str
    pre_velocity: Tuple[Vec3, Vec3, Vec3]
    mass_ratio: Tuple[float, float, float]
    expected_velocity: Optional[Tuple[Vec3, Vec3, Vec3]] = None
    body_to_field: Tuple[int, int, int] = (0, 1, 2)
    expected_contact: Optional[bool] = True
    friction: float = 0.0
    diagnostic: bool = False

    @property
    def case_id(self) -> str:
        return "threeFieldContact__" + self.key


@dataclass
class Tolerances:
    velocity: float = 1.0e-8
    momentum_abs: float = 1.0e-10
    momentum_rel: float = 1.0e-8
    impulse: float = 1.0e-10
    force: float = 1.0e-9
    mass: float = 1.0e-12
    gap: float = 1.0e-10
    canonical_velocity: float = 1.0e-5
    canonical_mass_relative: float = 0.10
    permutation: float = 1.0e-8


@dataclass
class Check:
    name: str
    status: str
    expected: str
    actual: str
    edge: str = ""
    required: bool = True


@dataclass
class Sample:
    state: int
    time: float
    mass: Tuple[float, float, float]
    active: Tuple[float, float, float]
    pre: Tuple[Vec3, Vec3, Vec3]
    post: Tuple[Vec3, Vec3, Vec3]
    force: Tuple[Vec3, Vec3, Vec3]
    position: Tuple[Vec3, Vec3, Vec3]
    normal: Tuple[Vec3, Vec3, Vec3]


@dataclass
class CaseResult:
    spec: CaseSpec
    status: str
    run_dir: str = ""
    history_path: str = ""
    selected_state: Optional[int] = None
    selected_time: Optional[float] = None
    checks: List[Check] = field(default_factory=list)
    samples: List[Sample] = field(default_factory=list, repr=False)
    plot_path: str = ""
    notes: List[str] = field(default_factory=list)


ZERO = (0.0, 0.0, 0.0)
SHARED_PRE = ((0.10, 0.0, 0.0), (0.05, 0.0, 0.0), (-0.10, 0.0, 0.0))
SHARED_POST = ((-0.0125, 0.0, 0.0),) * 3

CASE_SPECS: Tuple[CaseSpec, ...] = (
    CaseSpec(
        "threeField_collinearChain",
        "Collinear chain",
        "Baseline coupled constraints sharing the middle field",
        "chain",
        ((0.10, 0.0, 0.0), ZERO, (-0.10, 0.0, 0.0)),
        (0.5, 3.0, 0.5),
        (ZERO, ZERO, ZERO),
    ),
    CaseSpec(
        "threeField_collinearRedundant",
        "Collinear redundant constraints",
        "Redundant parallel A-C constraint and PGS ordering sensitivity",
        "chain",
        ((0.10, 0.0, 0.0), ZERO, (-0.10, 0.0, 0.0)),
        (0.5, 3.0, 0.5),
        (ZERO, ZERO, ZERO),
    ),
    CaseSpec(
        "threeField_collinearExplicitSurfaces",
        "Collinear explicit surfaces",
        "Known one-surface-position-per-field representation limitation",
        "chain_explicit",
        ((0.10, 0.0, 0.0), ZERO, (-0.10, 0.0, 0.0)),
        (0.5, 3.0, 0.5),
        (ZERO, ZERO, ZERO),
        diagnostic=True,
    ),
    CaseSpec(
        "threeField_sharedInterfaceOverlap",
        "Shared interface overlap",
        "Clean two-constraint coupled projection",
        "shared",
        SHARED_PRE,
        (1.0, 1.0, 2.0),
        SHARED_POST,
    ),
    CaseSpec(
        "threeField_sharedInterfaceZeroGap",
        "Shared interface zero gap",
        "Exact activation boundary at g = 0",
        "shared_zero",
        SHARED_PRE,
        (1.0, 1.0, 2.0),
        SHARED_POST,
    ),
    CaseSpec(
        "threeField_sharedInterfacePositiveGap",
        "Shared interface positive gap",
        "No premature contact for separated approaching bodies",
        "shared_no_contact",
        SHARED_PRE,
        (1.0, 1.0, 2.0),
        SHARED_PRE,
        expected_contact=False,
    ),
    CaseSpec(
        "threeField_sharedInterfaceSeparating",
        "Shared interface separating",
        "Unilateral complementarity for overlapping but separating bodies",
        "shared_no_contact",
        ((-0.10, 0.0, 0.0), (-0.05, 0.0, 0.0), (0.10, 0.0, 0.0)),
        (1.0, 1.0, 2.0),
        ((-0.10, 0.0, 0.0), (-0.05, 0.0, 0.0), (0.10, 0.0, 0.0)),
        expected_contact=False,
    ),
    CaseSpec(
        "threeField_sharedInterfacePermuted",
        "Shared interface permuted groups",
        "Contact-field numbering and constraint-order invariance",
        "shared",
        SHARED_PRE,
        (1.0, 1.0, 2.0),
        SHARED_POST,
        body_to_field=(2, 0, 1),
    ),
    CaseSpec(
        "threeField_sharedInterfaceMassRatio",
        "Shared interface mass ratio",
        "Ill-conditioned 1:100:2 nodal mass ratio",
        "shared",
        SHARED_PRE,
        (1.0, 100.0, 2.0),
        (((4.9 / 103.0), 0.0, 0.0),) * 3,
    ),
    CaseSpec(
        "threeField_sharedInterfaceFriction",
        "Shared interface friction",
        "Two Coulomb projections sharing the right-hand field",
        "shared_friction",
        ((0.10, 0.20, 0.0), (0.05, -0.20, 0.0), (-0.10, 0.0, 0.0)),
        (1.0, 1.0, 2.0),
        ((-0.0125, 0.16625, 0.0), (-0.0125, -0.18125, 0.0), (-0.0125, 0.0075, 0.0)),
        friction=0.30,
    ),
    CaseSpec(
        "threeField_corner",
        "Three-field corner",
        "Nonparallel competing normals at a triple junction",
        "corner",
        ((0.10, 0.10, 0.0), ZERO, ZERO),
        (1.0, 1.0, 1.0),
        None,
    ),
)
SPEC_BY_KEY = {spec.key: spec for spec in CASE_SPECS}
SPEC_BY_ID = {spec.case_id: spec for spec in CASE_SPECS}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--source-dir", default=str(Path(__file__).resolve().parent))
    parser.add_argument("--run-dir", default="", help="single-case harness mode")
    parser.add_argument("--runs-root", default="", help="fallback root containing case run directories")
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--case-id", default="")
    parser.add_argument("--case", action="append", default=[], help="case key; repeat to select cases")
    parser.add_argument("--output-prefix", default="threeFieldContact")
    parser.add_argument("--python", dest="python_cmd", default=sys.executable)
    parser.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", os.environ.get("VISIT_CMD", "")))
    parser.add_argument("--visit-timeout", type=float, default=300.0)
    parser.add_argument("--no-visit", action="store_true", help="do not extract missing CSV histories")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--force", "-y", action="store_true", help="accepted for harness compatibility")
    parser.add_argument("--job-id", default="", help="accepted for harness compatibility")
    parser.add_argument("--velocity-tol", type=float, default=1.0e-8)
    parser.add_argument("--momentum-abs-tol", type=float, default=1.0e-10)
    parser.add_argument("--momentum-rel-tol", type=float, default=1.0e-8)
    parser.add_argument("--impulse-tol", type=float, default=1.0e-10)
    parser.add_argument("--force-tol", type=float, default=1.0e-9)
    parser.add_argument("--gap-tol", type=float, default=1.0e-10)
    parser.add_argument("--permutation-tol", type=float, default=1.0e-8)
    return parser.parse_args(argv)


def finite(value: float) -> bool:
    return math.isfinite(value)


def vec_finite(vector: Vec3) -> bool:
    return all(finite(value) for value in vector)


def vsub(left: Vec3, right: Vec3) -> Vec3:
    return tuple(left[i] - right[i] for i in range(3))  # type: ignore[return-value]


def vnorm(vector: Vec3) -> float:
    return math.sqrt(sum(value * value for value in vector))


def vmax_abs(vectors: Iterable[Vec3]) -> float:
    values = [abs(value) for vector in vectors for value in vector if finite(value)]
    return max(values) if values else float("nan")


def format_number(value: float, digits: int = 6) -> str:
    if not finite(value):
        return "missing"
    return ("{0:." + str(digits) + "g}").format(value)


def parse_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def csv_value(row: Dict[str, str], name: str) -> float:
    aliases = (name, name.lower(), name.replace("grid", ""), name.replace("grid", "").lower())
    lowered = {key.lower(): value for key, value in row.items()}
    for alias in aliases:
        if alias in row:
            return parse_float(row[alias])
        if alias.lower() in lowered:
            return parse_float(lowered[alias.lower()])
    return float("nan")


def read_history(path: Path, spec: CaseSpec) -> List[Sample]:
    samples: List[Sample] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader):
            raw_mass = tuple(csv_value(row, "gridMass_f{0}".format(i)) for i in range(3))
            raw_active = tuple(csv_value(row, "gridActive_f{0}".format(i)) for i in range(3))

            def raw_vectors(base: str) -> Tuple[Vec3, Vec3, Vec3]:
                return tuple(
                    tuple(csv_value(row, "{0}_f{1}_{2}".format(base, field_index, component)) for component in COMPONENTS)
                    for field_index in range(3)
                )  # type: ignore[return-value]

            raw_pre = raw_vectors("gridUncontactedVelocity")
            raw_post = raw_vectors("gridVelocity")
            raw_force = raw_vectors("gridContactForce")
            raw_position = raw_vectors("gridSurfacePosition")
            raw_normal = raw_vectors("gridSurfaceNormal")
            mapping = spec.body_to_field
            parsed_state = parse_float(row.get("state", row_number))
            if not finite(parsed_state):
                parsed_state = float(row_number)
            samples.append(
                Sample(
                    state=int(parsed_state),
                    time=parse_float(row.get("time", row_number)),
                    mass=tuple(raw_mass[index] for index in mapping),
                    active=tuple(raw_active[index] for index in mapping),
                    pre=tuple(raw_pre[index] for index in mapping),
                    post=tuple(raw_post[index] for index in mapping),
                    force=tuple(raw_force[index] for index in mapping),
                    position=tuple(raw_position[index] for index in mapping),
                    normal=tuple(raw_normal[index] for index in mapping),
                )
            )
    return samples


def sample_has_core(sample: Sample, tolerances: Tolerances) -> bool:
    return all(finite(mass) and mass > tolerances.mass for mass in sample.mass) and all(vec_finite(vector) for vector in sample.post)


def has_pre(sample: Sample) -> bool:
    return all(vec_finite(vector) for vector in sample.pre)


def impulses(sample: Sample) -> Tuple[Vec3, Vec3, Vec3]:
    if not has_pre(sample):
        nan = (float("nan"),) * 3
        return (nan, nan, nan)  # type: ignore[return-value]
    return tuple(
        tuple(sample.mass[body] * (sample.post[body][component] - sample.pre[body][component]) for component in range(3))
        for body in range(3)
    )  # type: ignore[return-value]


def activity(sample: Sample) -> float:
    body_impulses = impulses(sample)
    if all(vec_finite(vector) for vector in body_impulses):
        return sum(vnorm(vector) for vector in body_impulses)
    values = [vnorm(vector) for vector in sample.force if vec_finite(vector)]
    return sum(values) if values else float("nan")


def momentum(sample: Sample, velocities: Sequence[Vec3]) -> Vec3:
    return tuple(sum(sample.mass[body] * velocities[body][component] for body in range(3)) for component in range(3))  # type: ignore[return-value]


def momentum_error(sample: Sample) -> float:
    if not has_pre(sample):
        return float("nan")
    return vnorm(vsub(momentum(sample, sample.post), momentum(sample, sample.pre)))


def normalized_mass_error(sample: Sample, expected: Sequence[float]) -> float:
    if not all(finite(value) and value > 0.0 for value in sample.mass):
        return float("inf")
    scale_actual = sum(sample.mass)
    scale_expected = sum(expected)
    actual = [value / scale_actual for value in sample.mass]
    target = [value / scale_expected for value in expected]
    return max(abs(actual[i] - target[i]) / max(target[i], 1.0e-30) for i in range(3))


def pre_velocity_error(sample: Sample, spec: CaseSpec) -> float:
    if not has_pre(sample):
        return float("inf")
    return max(abs(sample.pre[body][component] - spec.pre_velocity[body][component]) for body in range(3) for component in range(3))


def select_sample(samples: Sequence[Sample], spec: CaseSpec, tolerances: Tolerances) -> Optional[Sample]:
    eligible = [sample for sample in samples if sample_has_core(sample, tolerances)]
    if not eligible:
        return None
    canonical = sorted(
        eligible,
        key=lambda sample: (
            pre_velocity_error(sample, spec),
            normalized_mass_error(sample, spec.mass_ratio),
            sample.time if finite(sample.time) else sample.state,
        ),
    )
    if pre_velocity_error(canonical[0], spec) <= tolerances.canonical_velocity:
        return canonical[0]
    if spec.expected_contact:
        active_samples = [sample for sample in eligible if finite(activity(sample)) and activity(sample) > tolerances.impulse]
        if active_samples:
            return min(active_samples, key=lambda sample: sample.time if finite(sample.time) else sample.state)
    return min(eligible, key=lambda sample: sample.time if finite(sample.time) else sample.state)


def add_check(
    checks: List[Check],
    name: str,
    condition: Optional[bool],
    expected: str,
    actual: str,
    edge: str = "",
    required: bool = True,
) -> None:
    status = "SKIP" if condition is None else ("PASS" if condition else "FAIL")
    checks.append(Check(name, status, expected, actual, edge, required))


def scan_logs(run_dir: Path) -> Tuple[List[str], Optional[int], List[str]]:
    failure_pattern = re.compile(
        r"(projected\s*gauss.seidel.{0,100}(failed|did not converge|non.?converg)|"
        r"newton.raphson.{0,100}(failed|did not converge|non.?converg)|"
        r"contact.{0,80}(failed|did not converge)|segmentation fault|floating point (exception|error)|"
        r"mpi_abort|\bnan\b|\binf(?:inity)?\b|out_of_memory|time limit)",
        re.IGNORECASE,
    )
    iteration_pattern = re.compile(
        r"(?:CoupledContactSolverDiagnostics[^\n]*maximumNodeIterations=|"
        r"(?:PGS|ProjectedGaussSeidel|NewtonRaphson).{0,80}?iterations?\s*[:=]?\s*)(\d+)",
        re.IGNORECASE,
    )
    failures: List[str] = []
    iterations: List[int] = []
    aggregate_records = 0
    diagnostic_solvers = set()
    diagnostic_nonconverged: List[int] = []
    diagnostic_residuals: List[float] = []
    diagnostic_numerical_guards: List[int] = []
    diagnostic_skipped_nodes: List[int] = []
    diagnostic_skipped_pairs: List[int] = []
    diagnostic_newton_fallbacks: List[int] = []
    diagnostic_solver_rollbacks: List[int] = []
    detailed_nodes = 0
    nonconverged_detailed_nodes = 0
    pair_records = 0
    failure_nodes = 0
    failure_fields = 0
    failure_pairs = 0
    inactive_pairs = 0
    inactive_zero_gap_pairs = 0
    coulomb_margins: List[float] = []
    fitted_normals: List[str] = []

    def token(line: str, name: str) -> Optional[str]:
        match = re.search(r"(?:^|\s)" + re.escape(name) + r"=([^\s]+)", line)
        return match.group(1).rstrip(".") if match else None

    def integer_token(line: str, name: str) -> Optional[int]:
        value = token(line, name)
        try:
            return int(value) if value is not None else None
        except ValueError:
            return None

    def float_token(line: str, name: str) -> Optional[float]:
        value = token(line, name)
        try:
            parsed = float(value) if value is not None else None
        except ValueError:
            return None
        return parsed if parsed is not None and math.isfinite(parsed) else None

    paths = sorted(set(run_dir.glob("*.out")) | set(run_dir.glob("*.log")))
    for path in paths:
        if "postProcess" in path.name or "visit" in path.name.lower():
            continue
        try:
            text = path.read_text(errors="replace")[-2_000_000:]
        except OSError:
            continue
        for line in text.splitlines():
            if failure_pattern.search(line):
                failures.append(path.name + ": " + line.strip()[:300])
            if "CoupledContactSolverDiagnostics" in line:
                aggregate_records += 1
                solver = token(line, "solver")
                if solver:
                    diagnostic_solvers.add(solver)
                nonconverged = integer_token(line, "nonconvergedNodes")
                if nonconverged is not None:
                    diagnostic_nonconverged.append(nonconverged)
                    if nonconverged > 0:
                        failures.append(path.name + ": " + line.strip()[:300])
                residual = float_token(line, "maximumVelocityResidual")
                if residual is not None:
                    diagnostic_residuals.append(residual)
                for name, destination in (
                    ("numericalGuardActivations", diagnostic_numerical_guards),
                    ("numericallySkippedNodes", diagnostic_skipped_nodes),
                    ("numericallySkippedPairs", diagnostic_skipped_pairs),
                    ("newtonToPGSFallbackNodes", diagnostic_newton_fallbacks),
                    ("solverRollbackNodes", diagnostic_solver_rollbacks),
                ):
                    value = integer_token(line, name)
                    if value is not None:
                        destination.append(value)
            elif "CoupledContactFailureNodeDiagnostics" in line:
                failure_nodes += 1
                if integer_token(line, "converged") == 0:
                    nonconverged_detailed_nodes += 1
            elif "CoupledContactFailureFieldDiagnostics" in line:
                failure_fields += 1
            elif "CoupledContactFailurePairDiagnostics" in line:
                failure_pairs += 1
            elif "CoupledContactNodeDiagnostics" in line:
                detailed_nodes += 1
                if integer_token(line, "converged") == 0:
                    nonconverged_detailed_nodes += 1
            elif "CoupledContactPairDiagnostics" in line:
                pair_records += 1
                active = integer_token(line, "active")
                bilateral = integer_token(line, "bilateral")
                gap = float_token(line, "gap")
                if active == 0:
                    inactive_pairs += 1
                    if gap is not None and abs(gap) <= 1.0e-12:
                        inactive_zero_gap_pairs += 1
                if active == 1 and bilateral == 0:
                    margin = float_token(line, "coulombMargin")
                    if margin is not None:
                        coulomb_margins.append(margin)
                normal_match = re.search(r"(?:^|\s)normal=(\[[^\]]+\])", line)
                if normal_match and normal_match.group(1) not in fitted_normals and len(fitted_normals) < 4:
                    fitted_normals.append(normal_match.group(1))
        iterations.extend(int(match.group(1)) for match in iteration_pattern.finditer(text))

    diagnostic_notes: List[str] = []
    if aggregate_records:
        diagnostic_notes.append(
            "Coupled-solver log: solver={0}; records={1}; max node iterations={2}; "
            "max velocity residual={3}; max nonconverged nodes={4}.".format(
                ",".join(sorted(diagnostic_solvers)) or "unknown",
                aggregate_records,
                max(iterations) if iterations else "not reported",
                format_number(max(diagnostic_residuals)) if diagnostic_residuals else "not reported",
                max(diagnostic_nonconverged) if diagnostic_nonconverged else 0,
            )
        )
        diagnostic_notes.append(
            "Numerical safeguards: max guard activations={0}; skipped nodes={1}; "
            "skipped pairs={2}; Newton-to-PGS fallbacks={3}; solver rollbacks={4}.".format(
                max(diagnostic_numerical_guards) if diagnostic_numerical_guards else 0,
                max(diagnostic_skipped_nodes) if diagnostic_skipped_nodes else 0,
                max(diagnostic_skipped_pairs) if diagnostic_skipped_pairs else 0,
                max(diagnostic_newton_fallbacks) if diagnostic_newton_fallbacks else 0,
                max(diagnostic_solver_rollbacks) if diagnostic_solver_rollbacks else 0,
            )
        )
    if detailed_nodes or pair_records:
        diagnostic_notes.append(
            "Pair-frame log: detailed nodes={0} ({1} nonconverged); pairs={2}; "
            "inactive pairs={3}; inactive near-zero-gap pairs={4}; minimum active Coulomb margin={5}; "
            "sample fitted normals={6}.".format(
                detailed_nodes,
                nonconverged_detailed_nodes,
                pair_records,
                inactive_pairs,
                inactive_zero_gap_pairs,
                format_number(min(coulomb_margins)) if coulomb_margins else "not reported",
                ", ".join(fitted_normals) if fitted_normals else "not reported",
            )
        )
    if failure_nodes or failure_fields or failure_pairs:
        diagnostic_notes.append(
            "Failure-detail log: nodes={0}; fields={1}; pairs={2}.".format(
                failure_nodes,
                failure_fields,
                failure_pairs,
            )
        )
    return failures[:20], (max(iterations) if iterations else None), diagnostic_notes


def evaluate_case(spec: CaseSpec, samples: List[Sample], run_dir: Optional[Path], history_path: Optional[Path], tolerances: Tolerances) -> CaseResult:
    checks: List[Check] = []
    result = CaseResult(spec=spec, status="INCOMPLETE", run_dir=str(run_dir or ""), history_path=str(history_path or ""), samples=samples)
    if not samples:
        add_check(checks, "Simulation history available", False, "junction-node CSV with one or more states", "no history was found", spec.edge_case)
        result.checks = checks
        result.notes.append("Run the case and rerun this script with VisIt available to extract the Silo history.")
        return result

    core = [sample for sample in samples if sample_has_core(sample, tolerances)]
    add_check(checks, "Three active fields at junction", bool(core), "positive mass and finite velocity for A, B, C", "{0}/{1} usable states".format(len(core), len(samples)), spec.edge_case)
    if not core:
        result.checks = checks
        result.status = "FAIL"
        return result

    selected = select_sample(samples, spec, tolerances)
    assert selected is not None
    result.selected_state = selected.state
    result.selected_time = selected.time

    failures: List[str] = []
    max_iterations: Optional[int] = None
    log_diagnostic_notes: List[str] = []
    if run_dir and run_dir.is_dir():
        failures, max_iterations, log_diagnostic_notes = scan_logs(run_dir)
    result.notes.extend(log_diagnostic_notes)
    add_check(checks, "Solver completed without fatal marker", not failures, "no coupled-contact/non-finite/fatal marker", failures[0] if failures else "no fatal marker found", "Convergence and numerical stability")
    add_check(
        checks,
        "Coupled-solver iteration limit",
        None if max_iterations is None else max_iterations <= 200,
        "maximum reported iterations <= 200",
        "not reported in logs" if max_iterations is None else str(max_iterations),
        "Iteration cap",
        required=False,
    )
    if spec.key == "threeField_collinearChain":
        add_check(
            checks,
            "Coupled solve uses multiple iterations",
            None if max_iterations is None else max_iterations >= 2,
            "reported coupled-solver iteration count >= 2",
            "not reported in logs" if max_iterations is None else str(max_iterations),
            "Middle field participates in both constraints",
            required=max_iterations is not None,
        )

    pre_available = [sample for sample in core if has_pre(sample)]
    add_check(
        checks,
        "Pre-contact velocity captured",
        bool(pre_available),
        "gridUncontactedVelocity for all three fields",
        "{0}/{1} usable states".format(len(pre_available), len(core)),
        "Same-step impulse and conservation oracle",
    )

    finite_all = all(
        finite(value)
        for sample in core
        for vector in sample.post
        for value in vector
    )
    add_check(checks, "Finite projected velocities", finite_all, "all projected components finite", "finite" if finite_all else "NaN/Inf present", "Numerical robustness")

    if pre_available:
        errors = [momentum_error(sample) for sample in pre_available]
        worst_error = max(errors)
        scales = [max(vnorm(momentum(sample, sample.pre)), vnorm(momentum(sample, sample.post)), 1.0) for sample in pre_available]
        allowed = max(tolerances.momentum_abs + tolerances.momentum_rel * scale for scale in scales)
        add_check(
            checks,
            "Same-step momentum conservation",
            worst_error <= allowed,
            "||Σm v_post - Σm v_pre|| <= abs+rel tolerance",
            "worst {0}; allowed {1}".format(format_number(worst_error), format_number(allowed)),
            "Conservative contact impulse",
        )

    selected_impulses = impulses(selected)
    selected_activity = activity(selected)
    if spec.expected_contact is True and not spec.diagnostic:
        add_check(
            checks,
            "Contact response present",
            finite(selected_activity) and selected_activity > tolerances.impulse,
            "nonzero same-step contact impulse",
            "Σ|J| = {0}".format(format_number(selected_activity)),
            spec.edge_case,
        )
    elif spec.expected_contact is False:
        no_contact_samples = pre_available
        worst_impulse = max((activity(sample) for sample in no_contact_samples), default=float("nan"))
        worst_force = max((vnorm(vector) for sample in core for vector in sample.force if vec_finite(vector)), default=float("nan"))
        no_impulse = finite(worst_impulse) and worst_impulse <= tolerances.impulse
        no_force = (not finite(worst_force)) or worst_force <= tolerances.force
        add_check(
            checks,
            "No unilateral contact impulse",
            no_impulse and no_force,
            "J = 0 and contact force = 0 at every sampled state",
            "max Σ|J|={0}, max |F_c|={1}".format(format_number(worst_impulse), format_number(worst_force)),
            spec.edge_case,
        )
        worst_delta = max((vmax_abs(vsub(sample.post[i], sample.pre[i]) for i in range(3)) for sample in no_contact_samples), default=float("nan"))
        add_check(
            checks,
            "Projection leaves velocity unchanged",
            finite(worst_delta) and worst_delta <= tolerances.velocity,
            "v_post = v_pre",
            "max |Δv| = {0}".format(format_number(worst_delta)),
            spec.edge_case,
        )

    if spec.family in ("shared", "shared_zero", "shared_friction"):
        margin_ac = selected.post[2][0] - selected.post[0][0]
        margin_bc = selected.post[2][0] - selected.post[1][0]
        add_check(
            checks,
            "Projected normal inequalities",
            margin_ac >= -tolerances.velocity and margin_bc >= -tolerances.velocity,
            "vC_x-vA_x >= 0 and vC_x-vB_x >= 0",
            "A-C={0}, B-C={1}".format(format_number(margin_ac), format_number(margin_bc)),
            "Coupled nonpenetration",
        )
    elif spec.family == "chain":
        margin_ab = selected.post[1][0] - selected.post[0][0]
        margin_bc = selected.post[2][0] - selected.post[1][0]
        add_check(
            checks,
            "Projected chain inequalities",
            margin_ab >= -tolerances.velocity and margin_bc >= -tolerances.velocity,
            "vB_x-vA_x >= 0 and vC_x-vB_x >= 0",
            "A-B={0}, B-C={1}".format(format_number(margin_ab), format_number(margin_bc)),
            "Coupled/redundant constraints",
        )
    elif spec.family == "corner":
        margin_ab = selected.post[1][0] - selected.post[0][0]
        margin_ac = selected.post[2][1] - selected.post[0][1]
        add_check(
            checks,
            "Independent corner inequalities",
            margin_ab >= -tolerances.velocity and margin_ac >= -tolerances.velocity,
            "vB_x-vA_x >= 0 and vC_y-vA_y >= 0",
            "x(A-B)={0}, y(A-C)={1}".format(format_number(margin_ab), format_number(margin_ac)),
            "Nonparallel normals",
        )
        if all(vec_finite(vector) for vector in selected_impulses):
            response_x = abs(selected_impulses[0][0]) + abs(selected_impulses[1][0])
            response_y = abs(selected_impulses[0][1]) + abs(selected_impulses[2][1])
            add_check(
                checks,
                "Two-axis contact response",
                response_x > tolerances.impulse and response_y > tolerances.impulse,
                "nonzero x and y projected impulses",
                "x={0}, y={1}".format(format_number(response_x), format_number(response_y)),
                "Corner must not collapse to one averaged normal",
            )

    if spec.family in ("shared", "shared_zero", "chain") and has_pre(selected):
        total_mass = sum(selected.mass)
        com_x = sum(selected.mass[i] * selected.pre[i][0] for i in range(3)) / total_mass
        common_error = max(abs(selected.post[i][0] - com_x) for i in range(3))
        add_check(
            checks,
            "Coupled projection value",
            common_error <= tolerances.velocity,
            "all normal velocities equal measured center-of-mass value",
            "v_cm={0}, max error={1}".format(format_number(com_x), format_number(common_error)),
            "Exact local projection oracle",
        )

    if spec.family == "shared_friction" and has_pre(selected):
        cones = []
        for body in (0, 1):
            normal_impulse = abs(selected_impulses[body][0])
            tangent_impulse = abs(selected_impulses[body][1])
            cones.append((tangent_impulse, spec.friction * normal_impulse))
        worst_cone = max(tangent - limit for tangent, limit in cones)
        add_check(
            checks,
            "Coulomb cones",
            worst_cone <= tolerances.impulse + tolerances.velocity,
            "|J_t| <= μ|J_n| for A-C and B-C, μ=0.30",
            "; ".join("{0}<={1}".format(format_number(tangent), format_number(limit)) for tangent, limit in cones),
            "Coupled friction projection",
        )
        slip_pre = abs(selected.pre[0][1] - selected.pre[2][1]) + abs(selected.pre[1][1] - selected.pre[2][1])
        slip_post = abs(selected.post[0][1] - selected.post[2][1]) + abs(selected.post[1][1] - selected.post[2][1])
        add_check(
            checks,
            "Tangential slip reduced",
            slip_post <= slip_pre + tolerances.velocity,
            "summed tangential relative speed does not increase",
            "pre={0}, post={1}".format(format_number(slip_pre), format_number(slip_post)),
            "Dissipative friction",
        )

    canonical = (
        pre_velocity_error(selected, spec) <= tolerances.canonical_velocity
        and normalized_mass_error(selected, spec.mass_ratio) <= tolerances.canonical_mass_relative
    )
    if spec.expected_velocity is not None:
        expected_error = max(
            abs(selected.post[body][component] - spec.expected_velocity[body][component])
            for body in range(3)
            for component in range(3)
        )
        add_check(
            checks,
            "Published numeric oracle",
            expected_error <= tolerances.velocity if canonical and not spec.diagnostic else None,
            "golden junction velocity for prescribed initial state",
            (
                "not applied to diagnostic-only case"
                if spec.diagnostic
                else ("max error={0}".format(format_number(expected_error)) if canonical else "not applied: sampled pre-state/mass ratio is no longer canonical")
            ),
            "Expected-versus-actual value",
            required=canonical and not spec.diagnostic,
        )

    if spec.family == "shared_zero" and pre_available:
        canonical_states = [
            sample
            for sample in pre_available
            if pre_velocity_error(sample, spec) <= tolerances.canonical_velocity
            and normalized_mass_error(sample, spec.mass_ratio) <= tolerances.canonical_mass_relative
        ]
        earliest = min(canonical_states, key=lambda sample: sample.time if finite(sample.time) else sample.state) if canonical_states else None
        earliest_activity = activity(earliest) if earliest else float("nan")
        add_check(
            checks,
            "Zero-gap immediate activation",
            (finite(earliest_activity) and earliest_activity > tolerances.impulse) if earliest else None,
            "closing contact produces impulse in first captured canonical state",
            "first canonical-state Σ|J|={0}".format(format_number(earliest_activity)) if earliest else "no canonical initial state was captured",
            "Strict g<0 versus tolerance-aware g<=tol",
            required=earliest is not None,
        )

    if spec.family == "chain_explicit":
        gap_ab = selected.position[1][0] - selected.position[0][0]
        gap_bc = selected.position[2][0] - selected.position[1][0]
        missed_signature = (
            vec_finite(selected.position[0])
            and vec_finite(selected.position[1])
            and vec_finite(selected.position[2])
            and gap_ab > tolerances.gap
            and gap_bc > tolerances.gap
            and ((not finite(selected_activity)) or selected_activity <= tolerances.impulse)
        )
        ideal_projection = max(abs(selected.post[i][0]) for i in range(3)) <= tolerances.velocity
        add_check(
            checks,
            "Explicit-surface diagnostic classification",
            missed_signature or ideal_projection,
            "known positive-gap/no-impulse signature, or corrected ideal projection",
            "gAB={0}, gBC={1}, Σ|J|={2}; {3}".format(
                format_number(gap_ab),
                format_number(gap_bc),
                format_number(selected_activity),
                "known limitation reproduced" if missed_signature else ("ideal projection achieved" if ideal_projection else "unclassified response"),
            ),
            spec.edge_case,
            required=False,
        )

    result.checks = checks
    required_failures = [check for check in checks if check.required and check.status == "FAIL"]
    required_skips = [check for check in checks if check.required and check.status == "SKIP"]
    if required_failures:
        result.status = "FAIL"
    elif required_skips:
        result.status = "INCOMPLETE"
    elif spec.diagnostic:
        result.status = "DIAGNOSTIC"
    else:
        result.status = "PASS"
    return result


def locate_spec(value: str) -> Optional[CaseSpec]:
    if not value:
        return None
    if value in SPEC_BY_ID:
        return SPEC_BY_ID[value]
    short = value.split("__", 1)[-1]
    if short in SPEC_BY_KEY:
        return SPEC_BY_KEY[short]
    for spec in CASE_SPECS:
        if value in (spec.title, spec.key.replace("threeField_", "")):
            return spec
    return None


def read_manifest_run_dir(source_dir: Path, spec: CaseSpec) -> Optional[Path]:
    case_output = source_dir / "output" / spec.case_id
    manifests = sorted(case_output.glob("*_jobs.json"), key=lambda path: path.stat().st_mtime, reverse=True)
    for manifest in manifests:
        try:
            payload = json.loads(manifest.read_text())
            candidate = Path(payload["run_dir"]).expanduser()
        except (OSError, KeyError, ValueError, TypeError):
            continue
        if candidate.is_dir():
            return candidate.resolve()
    return None


def discover_run_dir(source_dir: Path, spec: CaseSpec, runs_root: Optional[Path]) -> Optional[Path]:
    manifest = read_manifest_run_dir(source_dir, spec)
    if manifest:
        return manifest
    candidates: List[Path] = []
    if runs_root:
        candidates.extend((runs_root / spec.case_id, runs_root / "verification" / spec.case_id, runs_root / spec.key))
    candidates.extend((source_dir / "output" / spec.case_id, source_dir / spec.case_id))
    for candidate in candidates:
        if not candidate.is_dir():
            continue
        if silo_exists(candidate) or any(candidate.glob("*_jobs.json")) or any(candidate.glob("*history.csv")):
            return candidate.resolve()
    return None


def discover_history(run_dir: Optional[Path], case_output: Path) -> Optional[Path]:
    candidates = [
        case_output / "three_field_contact_history.csv",
        case_output / "threeFieldContact_history.csv",
    ]
    if run_dir:
        candidates.extend((run_dir / "three_field_contact_history.csv", run_dir / "threeFieldContact_history.csv"))
    return next((candidate.resolve() for candidate in candidates if candidate.is_file()), None)


def detect_visit(requested: str) -> Optional[str]:
    candidates = (requested, os.environ.get("VISIT_COMMAND", ""), os.environ.get("VISIT_CMD", ""), "/usr/gapps/visit/bin/visit", "visit")
    for candidate in candidates:
        if not candidate:
            continue
        if os.path.isabs(candidate) and os.access(candidate, os.X_OK):
            return candidate
        found = shutil.which(candidate)
        if found:
            return found
    return None


def silo_exists(run_dir: Optional[Path]) -> bool:
    if not run_dir:
        return False
    silo = run_dir / "siloFiles"
    return silo.is_dir() and bool(list(silo.glob("mpm_cpdi_*")) or list(silo.glob("mpm_*")))


def extract_history(args: argparse.Namespace, source_dir: Path, run_dir: Path, output: Path) -> Tuple[Optional[Path], str]:
    if args.no_visit:
        return None, "VisIt extraction disabled by --no-visit"
    visit = detect_visit(args.visit_cmd)
    if not visit:
        return None, "VisIt executable not found"
    extractor = run_dir / "visitExtract_threeFieldContact.py"
    if not extractor.is_file():
        extractor = source_dir / "visitExtract_threeFieldContact.py"
    if not extractor.is_file():
        return None, "VisIt extractor script not found"
    if not silo_exists(run_dir):
        return None, "Silo database not found"
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [visit, "-nowin", "-cli", "-s", str(extractor), "--run-dir", str(run_dir), "--output", str(output)]
    try:
        process = subprocess.run(
            command,
            cwd=run_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=max(1.0, args.visit_timeout),
        )
    except subprocess.TimeoutExpired as error:
        return None, "VisIt extraction timed out after {0}s".format(args.visit_timeout)
    log_path = output.with_suffix(".visit.log")
    log_path.write_text(process.stdout + "\nreturncode={0}\n".format(process.returncode))
    if process.returncode != 0 or not output.is_file():
        return None, "VisIt extraction failed; see " + str(log_path)
    return output.resolve(), ""


def compare_permutation(results: Sequence[CaseResult], tolerances: Tolerances) -> None:
    by_key = {result.spec.key: result for result in results}
    baseline = by_key.get("threeField_sharedInterfaceOverlap")
    permuted = by_key.get("threeField_sharedInterfacePermuted")
    if not permuted:
        return
    if not baseline or baseline.selected_state is None or permuted.selected_state is None:
        add_check(
            permuted.checks,
            "Permutation equivalence to baseline",
            None,
            "remapped body histories agree with sharedInterfaceOverlap",
            "baseline or permuted history unavailable",
            "Field-order invariance",
        )
        if permuted.status == "PASS":
            permuted.status = "INCOMPLETE"
        return
    left = select_sample(baseline.samples, baseline.spec, tolerances)
    right = select_sample(permuted.samples, permuted.spec, tolerances)
    if left is None or right is None:
        return
    velocity_error = max(abs(left.post[i][c] - right.post[i][c]) for i in range(3) for c in range(3))
    left_impulse = impulses(left)
    right_impulse = impulses(right)
    impulse_error = max(abs(left_impulse[i][c] - right_impulse[i][c]) for i in range(3) for c in range(3))
    condition = velocity_error <= tolerances.permutation and impulse_error <= tolerances.permutation
    add_check(
        permuted.checks,
        "Permutation equivalence to baseline",
        condition,
        "remapped velocity and impulse differences <= {0}".format(format_number(tolerances.permutation)),
        "max Δv={0}, max ΔJ={1}".format(format_number(velocity_error), format_number(impulse_error)),
        "Field-order invariance",
    )
    if not condition:
        permuted.status = "FAIL"


def create_plot(result: CaseResult, figures_dir: Path) -> Optional[Path]:
    if not result.samples or result.selected_state is None:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:
        result.notes.append("Plot generation skipped: {0}".format(error))
        return None
    selected = select_sample(result.samples, result.spec, Tolerances())
    if selected is None:
        return None
    figures_dir.mkdir(parents=True, exist_ok=True)
    path = figures_dir / (result.spec.key + ".png")
    colors = ("#2d6cdf", "#ef8354", "#36a269")
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), constrained_layout=True)
    x = list(range(3))
    width = 0.25
    for component, axis in ((0, axes[0][0]), (1, axes[0][1])):
        pre_values = [selected.pre[i][component] for i in range(3)]
        post_values = [selected.post[i][component] for i in range(3)]
        axis.bar([value - width for value in x], pre_values, width, label="actual pre", color="#b8c7dc")
        axis.bar(x, post_values, width, label="actual post", color="#2d6cdf")
        if result.spec.expected_velocity is not None:
            expected = [result.spec.expected_velocity[i][component] for i in range(3)]
            axis.bar([value + width for value in x], expected, width, label="expected", color="#ef8354", alpha=0.8)
        axis.axhline(0.0, color="#333333", linewidth=0.8)
        axis.set_xticks(x, BODY_NAMES)
        axis.set_ylabel("velocity")
        axis.set_title(("Normal x" if component == 0 else "Tangential y") + " at selected state")
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize=8)

    usable = [sample for sample in result.samples if all(vec_finite(vector) for vector in sample.post)]
    times = [sample.time if finite(sample.time) else float(sample.state) for sample in usable]
    constraint_axis = axes[1][0]
    if result.spec.family == "corner":
        values_1 = [sample.post[1][0] - sample.post[0][0] for sample in usable]
        values_2 = [sample.post[2][1] - sample.post[0][1] for sample in usable]
        labels = ("vB_x-vA_x", "vC_y-vA_y")
    elif result.spec.family in ("chain", "chain_explicit"):
        values_1 = [sample.post[1][0] - sample.post[0][0] for sample in usable]
        values_2 = [sample.post[2][0] - sample.post[1][0] for sample in usable]
        labels = ("vB_x-vA_x", "vC_x-vB_x")
    else:
        values_1 = [sample.post[2][0] - sample.post[0][0] for sample in usable]
        values_2 = [sample.post[2][0] - sample.post[1][0] for sample in usable]
        labels = ("vC_x-vA_x", "vC_x-vB_x")
    constraint_axis.plot(times, values_1, marker="o", markersize=3, label=labels[0], color=colors[0])
    constraint_axis.plot(times, values_2, marker="s", markersize=3, label=labels[1], color=colors[1])
    constraint_axis.axhline(0.0, color="#333333", linewidth=0.8)
    constraint_axis.set_xlabel("time")
    constraint_axis.set_ylabel("relative velocity")
    if result.spec.expected_contact is False or result.spec.diagnostic:
        constraint_axis.set_title("Relative normal speed (inactive pairs may be < 0)")
    else:
        constraint_axis.set_title("Projected constraint margins (expected ≥ 0)")
    constraint_axis.grid(alpha=0.25)
    constraint_axis.legend(fontsize=8)

    diagnostic_axis = axes[1][1]
    activities = [activity(sample) for sample in usable]
    momentum_errors = [momentum_error(sample) for sample in usable]
    diagnostic_axis.plot(times, activities, marker="o", markersize=3, label="Σ|contact impulse|", color=colors[2])
    diagnostic_axis.plot(times, momentum_errors, marker="s", markersize=3, label="momentum error", color="#8f5da2")
    diagnostic_axis.axhline(0.0, color="#333333", linewidth=0.8)
    diagnostic_axis.set_xlabel("time")
    diagnostic_axis.set_ylabel("magnitude")
    diagnostic_axis.set_title("Contact activity and conservation")
    diagnostic_axis.grid(alpha=0.25)
    diagnostic_axis.legend(fontsize=8)

    fig.suptitle("{0} — {1} (state {2}, t={3})".format(result.spec.title, result.status, selected.state, format_number(selected.time)), fontsize=14)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path.resolve()


def result_dict(result: CaseResult) -> Dict[str, object]:
    return {
        "case": result.spec.key,
        "case_id": result.spec.case_id,
        "title": result.spec.title,
        "edge_case": result.spec.edge_case,
        "status": result.status,
        "run_dir": result.run_dir,
        "history_path": result.history_path,
        "selected_state": result.selected_state,
        "selected_time": result.selected_time,
        "plot_path": result.plot_path,
        "notes": result.notes,
        "checks": [asdict(check) for check in result.checks],
    }


def write_csv(results: Sequence[CaseResult], path: Path) -> None:
    with path.open("w", newline="") as handle:
        fieldnames = ("case", "case_status", "check", "check_status", "required", "edge_case", "expected", "actual")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            for check in result.checks:
                writer.writerow(
                    {
                        "case": result.spec.key,
                        "case_status": result.status,
                        "check": check.name,
                        "check_status": check.status,
                        "required": check.required,
                        "edge_case": check.edge,
                        "expected": check.expected,
                        "actual": check.actual,
                    }
                )


def relative_link(target: str, report_dir: Path) -> str:
    if not target:
        return ""
    try:
        return os.path.relpath(target, report_dir)
    except ValueError:
        return target


def write_markdown(results: Sequence[CaseResult], overall: str, path: Path) -> None:
    lines = [
        "# Three-field contact verification report",
        "",
        "**Overall status: {0}**".format(overall),
        "",
        "| Test | Edge case | Result | Selected state/time |",
        "|---|---|---:|---:|",
    ]
    for result in results:
        selected = "—" if result.selected_state is None else "{0} / {1}".format(result.selected_state, format_number(result.selected_time or 0.0))
        lines.append("| {0} | {1} | **{2}** | {3} |".format(result.spec.title, result.spec.edge_case, result.status, selected))
    for result in results:
        lines.extend(("", "## {0} — {1}".format(result.spec.title, result.status), "", result.spec.edge_case, ""))
        if result.plot_path:
            lines.extend(("![Expected-versus-actual diagnostics]({0})".format(relative_link(result.plot_path, path.parent)), ""))
        lines.extend(("| Check | Result | Expected | Actual |", "|---|---:|---|---|"))
        for check in result.checks:
            lines.append("| {0} | {1} | {2} | {3} |".format(check.name, check.status, check.expected.replace("|", "\\|"), check.actual.replace("|", "\\|")))
        for note in result.notes:
            lines.append("\n> " + note)
    path.write_text("\n".join(lines) + "\n")


def pdf_safe(value: object) -> str:
    """Return portable ASCII text while preserving mathematical meaning."""
    text = str(value)
    replacements = {
        "\u2011": "-",
        "\u2013": "-",
        "\u2014": "-",
        "\u2212": "-",
        "\u03a3": "Sum ",
        "\u03bc": "mu",
        "\u0394": "delta ",
        "\u2265": ">=",
        "\u2264": "<=",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text.encode("ascii", "replace").decode("ascii")


def write_pdf(results: Sequence[CaseResult], overall: str, path: Path) -> None:
    """Write the human-readable verification report as a paginated PDF.

    Matplotlib is already used for the expected-versus-actual plots, so using
    its PdfPages backend here avoids making ReportLab a separate runtime
    requirement for the postprocessor.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.lines import Line2D
        from matplotlib.patches import FancyBboxPatch
    except ImportError as error:
        raise RuntimeError("PDF report generation requires matplotlib") from error

    path.parent.mkdir(parents=True, exist_ok=True)
    status_colors = {
        "PASS": "#137A4B",
        "FAIL": "#B42318",
        "INCOMPLETE": "#9A6700",
        "DIAGNOSTIC": "#5B4BC4",
        "SKIP": "#667085",
    }
    ink = "#1E293B"
    muted = "#64748B"
    line = "#D8E0E8"
    pale = "#F5F8FB"
    white = "#FFFFFF"

    def wrapped(value: object, width: int) -> str:
        import textwrap

        text = pdf_safe(value)
        return "\n".join(textwrap.wrap(text, width=max(8, width), break_long_words=True, break_on_hyphens=False)) or " "

    def new_page():
        figure = plt.figure(figsize=(11.0, 8.5), facecolor=white)
        figure.subplots_adjust(0, 0, 1, 1)
        return figure

    def add_footer(figure, page_number: int) -> None:
        figure.add_artist(Line2D([0.055, 0.945], [0.052, 0.052], transform=figure.transFigure, color=line, linewidth=0.7))
        figure.text(0.055, 0.027, "GEOS MPM - Three-field contact verification", fontsize=6.7, color=muted, va="center")
        figure.text(0.945, 0.027, "Page {0}".format(page_number), fontsize=6.7, color=muted, va="center", ha="right")

    def style_table(table, statuses: Sequence[str], font_size: float, status_column: int) -> None:
        table.auto_set_font_size(False)
        table.set_fontsize(font_size)
        cells = table.get_celld()
        column_count = max(column for _, column in cells) + 1
        for column in range(column_count):
            cell = cells[(0, column)]
            cell.set_facecolor(ink)
            cell.set_edgecolor(line)
            cell.set_linewidth(0.6)
            cell.get_text().set_color(white)
            cell.get_text().set_weight("bold")
            cell.get_text().set_ha("left")
            cell.PAD = 0.035
        for row, status in enumerate(statuses, start=1):
            for column in range(column_count):
                cell = cells[(row, column)]
                cell.set_edgecolor(line)
                cell.set_linewidth(0.5)
                cell.set_facecolor(pale if row % 2 == 0 else white)
                cell.get_text().set_color(ink)
                cell.get_text().set_ha("left")
                cell.get_text().set_va("center")
                cell.PAD = 0.028
            status_cell = cells[(row, status_column)]
            status_cell.set_facecolor(status_colors.get(status, muted))
            status_cell.get_text().set_color(white)
            status_cell.get_text().set_weight("bold")
            status_cell.get_text().set_ha("center")

    counts = {status: sum(result.status == status for result in results) for status in ("PASS", "FAIL", "INCOMPLETE", "DIAGNOSTIC")}
    metadata = {
        "Title": "Three-field contact verification report",
        "Author": "GEOS MPM verification postprocessor",
        "Subject": "Expected-versus-actual projected contact verification",
    }
    with PdfPages(str(path), metadata=metadata) as pdf:
        page_number = 1
        figure = new_page()
        figure.text(0.055, 0.925, "Three-field contact verification report", fontsize=21, weight="bold", color=ink, va="top")
        figure.text(
            0.055,
            0.875,
            "Automated expected-versus-actual checks at the common grid node.",
            fontsize=9.2,
            color=muted,
            va="top",
        )
        figure.text(0.945, 0.922, "OVERALL  {0}".format(pdf_safe(overall)), fontsize=11, weight="bold", color=status_colors.get(overall, ink), ha="right", va="top")

        card_width = 0.205
        card_gap = 0.023
        for index, status in enumerate(("PASS", "FAIL", "INCOMPLETE", "DIAGNOSTIC")):
            x = 0.055 + index * (card_width + card_gap)
            card = FancyBboxPatch(
                (x, 0.755),
                card_width,
                0.082,
                boxstyle="round,pad=0.005,rounding_size=0.008",
                linewidth=0.7,
                edgecolor=line,
                facecolor=pale,
                transform=figure.transFigure,
            )
            figure.add_artist(card)
            figure.text(x + 0.014, 0.817, status, fontsize=7.2, weight="bold", color=muted, va="top")
            figure.text(x + 0.014, 0.785, str(counts[status]), fontsize=16, weight="bold", color=status_colors[status], va="center")

        figure.text(0.055, 0.716, "Suite summary", fontsize=13.5, weight="bold", color=ink, va="top")
        summary_rows = []
        summary_statuses = []
        for result in results:
            selected = "not available" if result.selected_state is None else "state {0}, t={1}".format(result.selected_state, format_number(result.selected_time or 0.0))
            summary_rows.append(
                [
                    wrapped(result.spec.title, 28),
                    wrapped(result.spec.edge_case, 68),
                    pdf_safe(result.status),
                    wrapped(selected, 25),
                ]
            )
            summary_statuses.append(result.status)
        summary_axis = figure.add_axes([0.055, 0.085, 0.89, 0.595])
        summary_axis.axis("off")
        summary_table = summary_axis.table(
            cellText=summary_rows,
            colLabels=["Test", "Edge case", "Result", "Selected sample"],
            colWidths=[0.19, 0.50, 0.11, 0.20],
            cellLoc="left",
            bbox=[0, 0, 1, 1],
        )
        style_table(summary_table, summary_statuses, 6.5, 2)
        add_footer(figure, page_number)
        pdf.savefig(figure, facecolor=white)
        plt.close(figure)

        for result in results:
            page_number += 1
            figure = new_page()
            figure.text(0.055, 0.930, pdf_safe(result.spec.title), fontsize=17, weight="bold", color=ink, va="top")
            figure.text(
                0.945,
                0.928,
                pdf_safe(result.status),
                fontsize=12,
                weight="bold",
                color=status_colors.get(result.status, muted),
                ha="right",
                va="top",
            )
            figure.text(0.055, 0.888, wrapped(result.spec.edge_case, 150), fontsize=8.4, color=ink, va="top")
            selected = "not available" if result.selected_state is None else "state {0}, t={1}".format(result.selected_state, format_number(result.selected_time or 0.0))
            history = wrapped(result.history_path or "not available", 145)
            figure.text(0.055, 0.833, "SELECTED SAMPLE", fontsize=6.4, weight="bold", color=muted, va="top")
            figure.text(0.170, 0.833, pdf_safe(selected), fontsize=7.2, color=ink, va="top")
            figure.text(0.055, 0.808, "HISTORY", fontsize=6.4, weight="bold", color=muted, va="top")
            figure.text(0.170, 0.808, history, fontsize=6.4, color=ink, va="top")

            has_plot = bool(result.plot_path and Path(result.plot_path).is_file())
            if has_plot:
                plot_axis = figure.add_axes([0.12, 0.435, 0.76, 0.345])
                plot_axis.imshow(plt.imread(str(result.plot_path)))
                plot_axis.axis("off")
                figure.text(
                    0.5,
                    0.422,
                    "Velocity bars show pre/post/oracle values; lower panels show constraint margins, contact activity, and momentum error.",
                    fontsize=6.3,
                    color=muted,
                    ha="center",
                    va="top",
                )
                table_height = min(0.305, max(0.160, 0.027 * (len(result.checks) + 1)))
                table_bottom = 0.385 - table_height
                table_font = 5.6 if len(result.checks) >= 9 else 6.1
            else:
                figure.text(0.055, 0.748, "No simulation figure is available for this case.", fontsize=7.2, color=muted, va="top")
                table_height = min(0.500, max(0.140, 0.055 * (len(result.checks) + 1)))
                table_bottom = 0.685 - table_height
                table_font = 6.7

            figure.text(0.055, table_bottom + table_height + 0.013, "Numerical checks", fontsize=9.0, weight="bold", color=ink, va="bottom")
            check_rows = []
            check_statuses = []
            for check in result.checks:
                check_rows.append(
                    [
                        wrapped(check.name, 28),
                        pdf_safe(check.status),
                        wrapped(check.expected, 58),
                        wrapped(check.actual, 63),
                    ]
                )
                check_statuses.append(check.status)
            if not check_rows:
                check_rows = [["No checks available", "INCOMPLETE", "Simulation output is required", "No history was found"]]
                check_statuses = ["INCOMPLETE"]
            check_axis = figure.add_axes([0.055, table_bottom, 0.89, table_height])
            check_axis.axis("off")
            check_table = check_axis.table(
                cellText=check_rows,
                colLabels=["Check", "Result", "Expected", "Actual"],
                colWidths=[0.18, 0.09, 0.35, 0.38],
                cellLoc="left",
                bbox=[0, 0, 1, 1],
            )
            style_table(check_table, check_statuses, table_font, 1)
            if result.notes:
                note_text = " | ".join("Note: " + pdf_safe(note) for note in result.notes)
                figure.text(0.055, 0.062, wrapped(note_text, 180), fontsize=5.8, color=muted, va="bottom")
            add_footer(figure, page_number)
            pdf.savefig(figure, facecolor=white)
            plt.close(figure)


def write_tex(results: Sequence[CaseResult], overall: str, path: Path) -> None:
    sanitize = lambda value: re.sub(r"[^A-Za-z0-9]", "", value)
    lines = ["% Generated by postProcess_threeFieldContact.py", "\\newcommand{\\threeFieldContactOverall}{%s}" % overall]
    for result in results:
        lines.append("\\newcommand{\\threeField%s}{%s}" % (sanitize(result.spec.key), result.status))
    path.write_text("\n".join(lines) + "\n")


def selected_specs(args: argparse.Namespace) -> List[CaseSpec]:
    if args.run_dir or args.case_id:
        spec = locate_spec(args.case_id or Path(args.run_dir).name)
        if not spec:
            raise SystemExit("Unknown three-field case: " + (args.case_id or Path(args.run_dir).name))
        return [spec]
    if not args.case:
        return list(CASE_SPECS)
    specs = []
    for value in args.case:
        spec = locate_spec(value)
        if not spec:
            raise SystemExit("Unknown case: " + value)
        specs.append(spec)
    return specs


def overall_status(results: Sequence[CaseResult]) -> str:
    if any(result.status == "FAIL" for result in results):
        return "FAIL"
    if any(result.status == "INCOMPLETE" for result in results):
        return "INCOMPLETE"
    return "PASS"


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    source_dir = Path(args.source_dir).expanduser().resolve()
    single_case = bool(args.run_dir)
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else (source_dir / "output" / "threeFieldContact_report").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    runs_root = Path(args.runs_root).expanduser().resolve() if args.runs_root else None
    tolerances = Tolerances(
        velocity=args.velocity_tol,
        momentum_abs=args.momentum_abs_tol,
        momentum_rel=args.momentum_rel_tol,
        impulse=args.impulse_tol,
        force=args.force_tol,
        gap=args.gap_tol,
        permutation=args.permutation_tol,
    )

    results: List[CaseResult] = []
    for spec in selected_specs(args):
        case_output = output_dir if single_case else source_dir / "output" / spec.case_id
        run_dir = Path(args.run_dir).expanduser().resolve() if single_case else discover_run_dir(source_dir, spec, runs_root)
        history = discover_history(run_dir, case_output)
        extraction_note = ""
        if history is None and run_dir:
            case_output.mkdir(parents=True, exist_ok=True)
            history, extraction_note = extract_history(args, source_dir, run_dir, case_output / "three_field_contact_history.csv")
        samples: List[Sample] = []
        if history:
            try:
                samples = read_history(history, spec)
            except Exception as error:
                extraction_note = "Could not read history: {0}".format(error)
        result = evaluate_case(spec, samples, run_dir, history, tolerances)
        if extraction_note:
            result.notes.append(extraction_note)
        results.append(result)

    if not single_case or len(results) > 1:
        compare_permutation(results, tolerances)
    if not args.no_plots:
        for result in results:
            plot = create_plot(result, figures_dir)
            if plot:
                result.plot_path = str(plot)

    overall = overall_status(results)
    payload = {
        "schema_version": 1,
        "overall_status": overall,
        "tolerances": asdict(tolerances),
        "cases": [result_dict(result) for result in results],
    }
    (output_dir / "three_field_contact_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_csv(results, output_dir / "three_field_contact_checks.csv")
    write_markdown(results, overall, output_dir / "three_field_contact_report.md")
    write_pdf(results, overall, output_dir / "three_field_contact_report.pdf")
    write_tex(results, overall, output_dir / "threeFieldContact_results.tex")

    print("Three-field contact status: " + overall)
    for result in results:
        print("  {0:<48} {1}".format(result.spec.key, result.status))
    print("Report: " + str(output_dir / "three_field_contact_report.pdf"))
    if overall == "FAIL":
        return 1
    if overall == "INCOMPLETE" and not args.allow_incomplete:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
