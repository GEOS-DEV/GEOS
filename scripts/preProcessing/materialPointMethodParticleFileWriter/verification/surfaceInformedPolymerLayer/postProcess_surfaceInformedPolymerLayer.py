#!/usr/bin/env python3
"""Post-process the surface-informed polymer layer verification.

The folder runs four variants: continuum/CZ and tension/compression.  The
quantitative comparison uses the same one-dimensional film models used by the
constitutive kernels: the continuum branch retains mean stress during the
J2 return, and the cohesive branch preserves the same volumetric normal stress
while returning only the deviatoric normal/shear part.  The VisIt render requests
show the loading-direction stress and equivalent plastic strain at four output
states for each variant when VisIt is available.
"""
from __future__ import annotations

import math
import os
from pathlib import Path
import sys

import numpy as np

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

#[pfw_dependency] pfw:pfw_materials.py
import pfw_materials as material_db

VARIANTS = [
    {
        "name": "continuum_tension",
        "label": "Continuum tension",
        "case_name": "surfaceInformedPolymerLayer_continuum_tension",
        "model_kind": "continuum",
        "final_global_strain": 0.10,
    },
    {
        "name": "continuum_compression",
        "label": "Continuum compression",
        "case_name": "surfaceInformedPolymerLayer_continuum_compression",
        "model_kind": "continuum",
        "final_global_strain": -0.08,
    },
    {
        "name": "cohesive_tension",
        "label": "Cohesive-zone tension",
        "case_name": "surfaceInformedPolymerLayer_cohesive_tension",
        "model_kind": "cohesive",
        "final_global_strain": 0.10,
    },
    {
        "name": "cohesive_compression",
        "label": "Cohesive-zone compression",
        "case_name": "surfaceInformedPolymerLayer_cohesive_compression",
        "model_kind": "cohesive",
        "final_global_strain": -0.08,
    },
]
DEFAULT_NAME = "surfaceInformedPolymerLayer"

# Material and geometry constants from pfw_input_surfaceInformedPolymerLayer.py.
# The material database uses the MPM validation unit system (mm, us, mg, K),
# where stresses supplied to GEOS material models are in GPa.  The post-processor
# reports stresses in MPa.  It also reconstructs the x-axis from the actual
# deformation-gradient history written by GEOS when available, rather than from
# wall-clock time.  This matters because the verification input uses a cosine
# loading table and because scheduler restarts can leave nonuniform history
# spacing.
# Must match the slow cosine ramp in pfw_input_surfaceInformedPolymerLayer.py.
# It is used only when a history row does not include F11/length_y.
STOP_TIME = 20.0
COHESIVE_GAGE_LENGTH = 1.0
FILM_THICKNESS = 0.10
CONTINUUM_GAGE_LENGTH = 1.0
DOMAIN_WIDTH = FILM_THICKNESS
DOMAIN_DEPTH = 0.01
DOMAIN_VOLUME = DOMAIN_WIDTH * CONTINUUM_GAGE_LENGTH * DOMAIN_DEPTH
# Fallback area for reaction histories when an older run does not write length_x/length_z.
# Current inputs write both lengths, so this is only a defensive default.
CROSS_SECTION_AREA = DOMAIN_WIDTH * 0.01
STRESS_TO_MPA = 1000.0
# Analytical comparison properties are read from the same generic FKM/Viton validation
# material card used by pfw_input_surfaceInformedPolymerLayer.py.  Keep material data in
# pfw_materials.py rather than embedding a separate calibration data set in this
# post-processor.
ANALYTICAL_MATERIAL = material_db.vitonFKM75SurfacePolymer.copy()

POLYMER_BULK_MODULUS = float(ANALYTICAL_MATERIAL["defaultBulkModulus"])
POLYMER_SHEAR_MODULUS = float(ANALYTICAL_MATERIAL["defaultShearModulus"])
POLYMER_CONSTRAINED_MODULUS = POLYMER_BULK_MODULUS + 4.0 * POLYMER_SHEAR_MODULUS / 3.0
YIELD_STRENGTH = float(ANALYTICAL_MATERIAL["defaultYieldStrength"])
SOFTENING_MAGNITUDE = float(ANALYTICAL_MATERIAL["shearSofteningMagnitude"])
SOFTENING_R1 = float(ANALYTICAL_MATERIAL["shearSofteningShapeParameter1"])
SOFTENING_R2 = float(ANALYTICAL_MATERIAL["shearSofteningShapeParameter2"])
HARDENING_SLOPE = float(ANALYTICAL_MATERIAL["strainHardeningSlope"])
PRESSURE_ASYMMETRY = float(ANALYTICAL_MATERIAL.get("pressureAsymmetryAmplitude", 0.0))
MAXIMUM_STRETCH = float(ANALYTICAL_MATERIAL["maximumStretch"])


def final_global_strain(variant_name: str) -> float:
    variant = next((v for v in VARIANTS if v["name"] == variant_name), None)
    return float(variant.get("final_global_strain", 0.0)) if variant else 0.0


def gage_length(model_kind: str) -> float:
    """Return the distance over which the fTable strain is applied."""
    return CONTINUUM_GAGE_LENGTH if model_kind == "continuum" else COHESIVE_GAGE_LENGTH


def imposed_film_strain(model_kind: str, global_strain: float) -> float:
    """Reconstruct the local polymer film strain from imposed gage strain.

    The continuum verification now loads a finite-thickness polymer layer through
    elastic bars, while the cohesive-zone verification uses elastic blocks around a
    finite-thickness CZ.  In both branches the fTable strain is applied over the
    one-unit bar/block gage length and converted to a local film strain.
    """
    return global_strain * gage_length(model_kind) / FILM_THICKNESS



def softening(kappa: float) -> float:
    ratio = max(kappa, 0.0) / max(SOFTENING_R1, 1.0e-30)
    return SOFTENING_MAGNITUDE * math.exp(-ratio ** SOFTENING_R2)


def chain_stretch_from_film_state(film_strain: float, tangential_strain: float = 0.0) -> float:
    stretch = max(1.0 + film_strain, 1.0e-16)
    F = np.array([[1.0, tangential_strain, 0.0], [0.0, stretch, 0.0], [0.0, 0.0, 1.0]], dtype=float)
    principal_stretch = float(math.sqrt(max(np.linalg.eigvalsh(F.T @ F).max(), 0.0)))
    J = max(float(np.linalg.det(F)), 1.0e-16)
    return max(principal_stretch, principal_stretch / (J ** (1.0 / 3.0)))


def stretch_hardening(film_strain: float, tangential_strain: float = 0.0) -> float:
    lam = min(chain_stretch_from_film_state(film_strain, tangential_strain), MAXIMUM_STRETCH)
    lam = max(lam, 1.0)
    return HARDENING_SLOPE * (lam * lam - 1.0 / lam)


def stress_decomposition(stress: list[float]) -> tuple[float, float, list[float]]:
    """Return mean stress, von-Mises stress, and normalized deviator."""
    p = (stress[0] + stress[1] + stress[2]) / 3.0
    dev = [stress[0] - p, stress[1] - p, stress[2] - p]
    norm = math.sqrt(max(sum(v * v for v in dev), 0.0))
    q = math.sqrt(1.5) * norm
    if norm <= 1.0e-30:
        return p, 0.0, [0.0, 0.0, 0.0]
    return p, q, [v / norm for v in dev]


def stress_recomposition(p: float, q: float, direction: list[float]) -> list[float]:
    scale = math.sqrt(2.0 / 3.0) * q
    return [p + scale * direction[i] for i in range(3)]


def strain_magnitude(strain: list[float]) -> float:
    return math.sqrt(max(sum(v * v for v in strain), 0.0))


def elastic_increment_uniaxial_strain(delta_film_strain: float) -> list[float]:
    """Elastic stress increment for eps=[0, delta_eps_y, 0]."""
    K = POLYMER_BULK_MODULUS
    G = POLYMER_SHEAR_MODULUS
    e = delta_film_strain
    return [
        (K - 2.0 * G / 3.0) * e,
        (K + 4.0 * G / 3.0) * e,
        (K - 2.0 * G / 3.0) * e,
    ]


def elastic_strain_increment_from_stress_increment(delta_stress: list[float]) -> list[float]:
    """Isotropic compliance for normal components only."""
    K = POLYMER_BULK_MODULUS
    G = POLYMER_SHEAR_MODULUS
    dp, dq, direction = stress_decomposition(delta_stress)
    strain = [dp / (3.0 * K) for _ in range(3)]
    if dq > 0.0:
        # This is the normal-component equivalent of the C++ helper:
        # delta_eps_dev = sqrt(2/3) * delta_q * n / (2G).
        factor = math.sqrt(2.0 / 3.0) * dq / (2.0 * G)
        for i in range(3):
            strain[i] += factor * direction[i]
    return strain


def continuum_update(total_film_strain: float,
                     previous_film_strain: float,
                     stress_old: list[float],
                     plastic_strain_old: list[float],
                     kappa_old: float,
                     hardening_film_strain: float | None = None) -> tuple[list[float], list[float], float]:
    """One-dimensional continuum update for the uniaxial-strain film.

    The continuum kernel is a J2 radial return with retained mean stress.  The
    analytical comparison therefore cannot simply limit sigma_yy by the flow
    stress.  It must preserve the hydrostatic part of the constrained film
    response and return only the deviatoric invariant q to the flow surface.
    """
    delta_film_strain = total_film_strain - previous_film_strain
    trial = [stress_old[i] + elastic_increment_uniaxial_strain(delta_film_strain)[i] for i in range(3)]
    p_trial, q_trial, direction = stress_decomposition(trial)
    plastic_increment = [0.0, 0.0, 0.0]
    hardening_strain = total_film_strain if hardening_film_strain is None else hardening_film_strain
    flow = YIELD_STRENGTH + softening(kappa_old) + stretch_hardening(hardening_strain)
    stress_new = trial[:]
    plastic_new = plastic_strain_old[:]
    kappa_new = kappa_old
    for _ in range(32):
        yield_measure = q_trial + PRESSURE_ASYMMETRY * p_trial
        if yield_measure <= flow or q_trial <= 1.0e-16:
            stress_new = trial[:]
            plastic_new = plastic_strain_old[:]
            kappa_new = kappa_old
            break
        returned_q = min(q_trial, max(flow - PRESSURE_ASYMMETRY * p_trial, 0.0))
        stress_candidate = stress_recomposition(p_trial, returned_q, direction)
        delta_stress = [stress_candidate[i] - stress_old[i] for i in range(3)]
        elastic_strain_increment = elastic_strain_increment_from_stress_increment(delta_stress)
        total_strain_increment = [0.0, delta_film_strain, 0.0]
        plastic_increment = [total_strain_increment[i] - elastic_strain_increment[i] for i in range(3)]
        plastic_candidate = [plastic_strain_old[i] + plastic_increment[i] for i in range(3)]
        kappa_candidate = strain_magnitude(plastic_candidate)
        flow_candidate = YIELD_STRENGTH + softening(kappa_candidate) + stretch_hardening(hardening_strain)
        stress_new = stress_candidate
        plastic_new = plastic_candidate
        kappa_new = kappa_candidate
        if abs(flow_candidate - flow) <= 1.0e-10 * max(1.0, abs(flow_candidate)):
            break
        flow = flow_candidate
    return stress_new, plastic_new, kappa_new


def cohesive_update(film_strain: float,
                    plastic_normal_old: float,
                    plastic_tangential_old: float,
                    kappa_old: float,
                    tangential_strain: float = 0.0) -> tuple[float, float, float, float]:
    """Reduced normal-shear cohesive update matching the C++ thin-film law.

    The cohesive-zone projection now uses the same volumetric/deviatoric split as
    the continuum uniaxial-strain update.  The mean normal film stress p=K eps_n
    is retained, and only the deviatoric normal stress plus shear stress are
    returned to the scalar flow surface.  This prevents the analytical cohesive
    curve from being plotted on a different stress scale than the continuum film
    when the polymer is nearly incompressible.
    """
    p_trial = POLYMER_BULK_MODULUS * film_strain
    normal_deviatoric_modulus = 4.0 * POLYMER_SHEAR_MODULUS / 3.0
    s_n_trial = normal_deviatoric_modulus * (film_strain - plastic_normal_old)
    tau_trial = POLYMER_SHEAR_MODULUS * (tangential_strain - plastic_tangential_old)
    q_trial = math.sqrt(max(2.25 * s_n_trial * s_n_trial + 3.0 * tau_trial * tau_trial, 0.0))

    kappa = kappa_old
    s_n = s_n_trial
    tau = tau_trial
    sigma_n = p_trial + s_n
    plastic = False
    for _ in range(16):
        flow = YIELD_STRENGTH + softening(kappa) + stretch_hardening(film_strain, tangential_strain)
        yield_measure = q_trial + PRESSURE_ASYMMETRY * p_trial
        if yield_measure <= flow or q_trial <= 1.0e-16:
            return sigma_n, plastic_normal_old, plastic_tangential_old, kappa_old
        plastic = True
        returned_q = min(q_trial, max(flow - PRESSURE_ASYMMETRY * p_trial, 0.0))
        scale = returned_q / max(q_trial, 1.0e-16)
        s_n_new = scale * s_n_trial
        tau_new = scale * tau_trial
        delta_kappa = max(0.0, (q_trial - returned_q) / (3.0 * max(POLYMER_SHEAR_MODULUS, 1.0e-16)))
        kappa_new = kappa_old + delta_kappa
        s_n = s_n_new
        tau = tau_new
        sigma_n = p_trial + s_n
        if abs(kappa_new - kappa) <= 1.0e-10 * max(1.0, abs(kappa_new)):
            kappa = kappa_new
            break
        kappa = kappa_new
    if not plastic:
        return sigma_n, plastic_normal_old, plastic_tangential_old, kappa_old
    plastic_normal_new = film_strain - s_n / max(normal_deviatoric_modulus, 1.0e-16)
    plastic_tangential_new = tangential_strain - tau / max(POLYMER_SHEAR_MODULUS, 1.0e-16)
    return sigma_n, plastic_normal_new, plastic_tangential_new, max(kappa_old, kappa)

def geos_card_analytical_series(model_kind: str, final_strain: float, n: int = 251) -> tuple[list[float], list[float], list[float]]:
    xs: list[float] = []
    stresses_mpa: list[float] = []
    kappas: list[float] = []
    previous_film_strain = 0.0
    if model_kind == "continuum":
        stress = [0.0, 0.0, 0.0]
        plastic_strain = [0.0, 0.0, 0.0]
        kappa = 0.0
        previous_log_film_strain = 0.0
        for i in range(n):
            tfrac = i / max(n - 1, 1)
            film_strain = imposed_film_strain(model_kind, final_strain * tfrac)
            film_stretch = max(1.0 + film_strain, 1.0e-16)
            log_film_strain = math.log(film_stretch)
            stress, plastic_strain, kappa = continuum_update(
                log_film_strain,
                previous_log_film_strain,
                stress,
                plastic_strain,
                kappa,
                hardening_film_strain=film_strain,
            )
            previous_log_film_strain = log_film_strain
            xs.append(film_strain)
            stresses_mpa.append(stress[1] * STRESS_TO_MPA)
            kappas.append(kappa)
    else:
        plastic_normal = 0.0
        plastic_tangential = 0.0
        kappa = 0.0
        for i in range(n):
            tfrac = i / max(n - 1, 1)
            film_strain = imposed_film_strain(model_kind, final_strain * tfrac)
            stress_scalar, plastic_normal, plastic_tangential, kappa = cohesive_update(film_strain, plastic_normal, plastic_tangential, kappa)
            xs.append(film_strain)
            stresses_mpa.append(stress_scalar * STRESS_TO_MPA)
            kappas.append(kappa)
    return xs, stresses_mpa, kappas


def analytical_series(model_kind: str, final_strain: float, n: int = 251) -> tuple[list[float], list[float], list[float]]:
    return geos_card_analytical_series(model_kind, final_strain, n)


def pick_stress_y(row: dict) -> float | None:
    labels = ["yy", "22", "_1", "y"]
    for key in row:
        low = key.lower()
        if "stress" in low and any(tok in low for tok in labels):
            try:
                return float(row[key])
            except Exception:
                pass
    for key in row:
        if "stress" in key.lower():
            try:
                return float(row[key])
            except Exception:
                pass
    return None


def pick_plastic(row: dict) -> float | None:
    for key in row:
        low = key.lower()
        if "plastic" in low and any(tok in low for tok in ["strain", "eps", "magnitude", "eq"]):
            try:
                return float(row[key])
            except Exception:
                pass
    clean = {str(k).strip().lower().replace(" ", ""): v for k, v in row.items()}
    component_keys = ["epxx", "epyy", "epzz", "epyz", "epxz", "epxy"]
    components = []
    found_component = False
    for key in component_keys:
        try:
            components.append(float(clean[key]))
            found_component = True
        except Exception:
            components.append(0.0)
    if found_component:
        return math.sqrt(sum(value * value for value in components))
    return None


def box_material_volume(row: dict) -> float | None:
    for key, value in row.items():
        if str(key).strip().lower() == "volume":
            try:
                volume = float(value)
                return volume if volume > 0.0 else None
            except Exception:
                return None
    return None


def pick_box_material_stress_y(row: dict) -> float | None:
    stress = pick_stress_y(row)
    volume = box_material_volume(row)
    if stress is None:
        return None
    if volume is None:
        return stress
    # Existing GEOS box-average stress rows are stress*materialVolume/domainVolume.
    # Correct them here so reports made from old runs still show the material average.
    return stress * DOMAIN_VOLUME / max(volume, 1.0e-30)


def history_time(row: dict) -> float:
    try:
        return float(row.get("time", row.get("Time", row.get("t", 0.0))))
    except Exception:
        return 0.0


def cosine_load_fraction(time: float) -> float:
    """Cosine interpolation fraction used by the verification fTable."""
    x = max(0.0, min(time / STOP_TIME, 1.0))
    return 0.5 - 0.5 * math.cos(math.pi * x)


def global_strain_from_history(model_kind: str, row: dict, final_strain: float, time: float) -> float:
    """Return the imposed global y strain from a history row.

    reactionHistory.csv contains F11 in current GEOS-MPM builds.  Some older or
    partially generated runs may only contain length_y.  Falling back to the
    cosine load table keeps the post-processor usable for dry runs and incomplete
    histories, but the preferred x-coordinate is always the actual deformation
    written by the solver.
    """
    clean = {str(k).strip(): v for k, v in row.items()}
    for key in ["F11", "f11", "F_11", "Fyy", "FYY"]:
        if key in clean:
            try:
                return float(clean[key]) - 1.0
            except Exception:
                pass
    for key in ["length_y", "lengthY", "Ly"]:
        if key in clean:
            try:
                return float(clean[key]) / max(gage_length(model_kind), 1.0e-30) - 1.0
            except Exception:
                pass
    return final_strain * cosine_load_fraction(time)


def film_strain_from_history(model_kind: str, row: dict, final_strain: float) -> tuple[float, float]:
    time = history_time(row)
    return time, imposed_film_strain(model_kind, global_strain_from_history(model_kind, row, final_strain, time))


def domain_stretch_from_history(model_kind: str, row: dict, final_strain: float) -> tuple[float, float]:
    time = history_time(row)
    return time, 1.0 + global_strain_from_history(model_kind, row, final_strain, time)


def reaction_stress_from_row(row: dict) -> float | None:
    clean = {str(k).strip(): v for k, v in row.items()}
    ry_minus = clean.get("Ry-", clean.get("reactionYMinus", None))
    ry_plus = clean.get("Ry+", clean.get("reactionYPlus", None))
    values = []
    for value in [ry_minus, ry_plus]:
        try:
            values.append(float(value))
        except Exception:
            pass
    if not values:
        return None
    # Boundary reactions have opposite signs on the two y faces.  Use the
    # sign of the top reaction when possible so tension is positive and
    # compression is negative in the verification plot.
    signed_force = values[-1] if len(values) > 1 else values[0]
    if len(values) > 1 and abs(values[0]) > abs(values[1]):
        signed_force = -values[0]
    try:
        area = abs(float(clean.get("length_x", 0.0)) * float(clean.get("length_z", 0.0)))
    except Exception:
        area = 0.0
    if area <= 1.0e-30:
        area = CROSS_SECTION_AREA
    return signed_force / max(area, 1.0e-30)


def read_reaction_stress_y(run_dir: Path, model_kind: str, final_strain: float) -> list[tuple[float, float, float]]:
    path = find_first(run_dir, ["reactionHistory.csv"])
    if path is None:
        return []
    rows = read_csv_numeric(path)
    out: list[tuple[float, float, float]] = []
    for row in rows:
        clean = {str(k).strip(): v for k, v in row.items()}
        time, film_strain = film_strain_from_history(model_kind, clean, final_strain)
        stress = reaction_stress_from_row(clean)
        if stress is not None:
            out.append((time, film_strain, stress))
    out.sort(key=lambda item: item[0])
    filtered: list[tuple[float, float, float]] = []
    last_t = -math.inf
    for time, film_strain, stress in out:
        if time > last_t:
            filtered.append((time, film_strain, stress))
            last_t = time
    return filtered


def read_reaction_stress_domain_stretch_y(run_dir: Path, model_kind: str, final_strain: float) -> list[tuple[float, float, float]]:
    path = find_first(run_dir, ["reactionHistory.csv"])
    if path is None:
        return []
    rows = read_csv_numeric(path)
    out: list[tuple[float, float, float]] = []
    for row in rows:
        clean = {str(k).strip(): v for k, v in row.items()}
        time, domain_stretch = domain_stretch_from_history(model_kind, clean, final_strain)
        stress = reaction_stress_from_row(clean)
        if stress is not None:
            out.append((time, domain_stretch, stress))
    out.sort(key=lambda item: item[0])
    filtered: list[tuple[float, float, float]] = []
    last_t = -math.inf
    for time, domain_stretch, stress in out:
        if time > last_t:
            filtered.append((time, domain_stretch, stress))
            last_t = time
    return filtered


def interpolate_series(xs: list[float], ys: list[float], x: float) -> float:
    if not xs:
        return 0.0
    pairs = sorted(zip(xs, ys), key=lambda item: item[0])
    if x <= pairs[0][0]:
        return pairs[0][1]
    if x >= pairs[-1][0]:
        return pairs[-1][1]
    for (x0, y0), (x1, y1) in zip(pairs[:-1], pairs[1:]):
        if x0 <= x <= x1:
            if abs(x1 - x0) <= 1.0e-30:
                return y1
            a = (x - x0) / (x1 - x0)
            return (1.0 - a) * y0 + a * y1
    return pairs[-1][1]


def expected_at_film_strain(model_kind: str, final_strain: float, film_strain: float) -> tuple[float, float, float]:
    xs, ss, ps = analytical_series(model_kind, final_strain, n=501)
    return film_strain, interpolate_series(xs, ss, film_strain), interpolate_series(xs, ps, film_strain)


def expected_at_time(model_kind: str, final_strain: float, time: float) -> tuple[float, float, float]:
    film_strain = imposed_film_strain(model_kind, final_strain * cosine_load_fraction(time))
    return expected_at_film_strain(model_kind, final_strain, film_strain)


def normalize_measured_stress_units(measured: list[float], expected_mpa: list[float]) -> tuple[list[float], float]:
    """Return stresses in MPa and the scale applied to the measured series.

    GEOS writes material stresses and force/area reaction reductions in the MPM
    stress unit, GPa, for this verification.  Keep the conversion explicit
    instead of inferring units from the expected curve magnitude; the inference
    hid large cohesive-zone reaction errors when the wrong reaction measure was
    being plotted.
    """
    return [STRESS_TO_MPA * value for value in measured], STRESS_TO_MPA


def series_completion_status(model_kind: str, final_strain: float, measured_x: list[float]) -> tuple[str, bool, float | None, float]:
    target = imposed_film_strain(model_kind, final_strain)
    if not measured_x:
        return "missing", False, None, target
    final_x = measured_x[-1]
    tolerance = max(1.0e-3, 2.0e-2 * abs(target))
    complete = abs(final_x - target) <= tolerance
    return ("complete" if complete else "incomplete"), complete, final_x, target


def reaction_measure(run_dir: Path) -> str:
    xml_path = find_first(run_dir, ["mpm_surfaceInformedPolymerLayer.xml"])
    if xml_path is None:
        return "unknown"
    try:
        text = xml_path.read_text(errors="ignore")
    except Exception:
        return "unknown"
    if "useInternalForceAsFaceReaction=\"1\"" in text:
        return "internal_force"
    if "reactionHistory=\"1\"" in text:
        return "velocity_correction"
    return "unknown"


def compact_or_na(value, precision: int = 6) -> str:
    text = compact_float(value, precision)
    return text if text else "--"


def latex_plot_block(output_dir: Path, name: str, caption: str) -> list[str]:
    if not (output_dir / name).is_file():
        return []
    return [
        r"\begin{center}",
        r"\includegraphics[width=0.92\linewidth]{\CaseOutputDir/" + name + "}",
        r"{\footnotesize " + latex_escape(caption) + r"\par}",
        r"\end{center}",
    ]


def plot_metric(output_dir: Path, file_name: str, title: str, x_label: str, y_label: str, series: list[tuple[str, list[float], list[float]]]) -> str | None:
    series = [(label, x, y) for label, x, y in series if x and y]
    if not series:
        return None

    def draw_order(item: tuple[str, list[float], list[float]]) -> int:
        label = item[0].lower()
        is_analytical = label.startswith("analytical")
        is_cohesive = "cohesive" in label
        return (0 if is_analytical else 2) + (1 if is_cohesive else 0)

    def series_color(label: str) -> str | None:
        lower = label.lower()
        if "cohesive" in lower:
            return "tab:cyan" if "tension" in lower else "tab:red"
        if "continuum" in lower:
            return "tab:blue" if "tension" in lower else "tab:orange"
        return None

    plt = matplotlib()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for label, x, y in sorted(series, key=draw_order):
        lower = label.lower()
        is_analytical = lower.startswith("analytical")
        is_cohesive = "cohesive" in lower
        ax.plot(
            x,
            y,
            label=label,
            color=series_color(label),
            linestyle="--" if is_cohesive else "-",
            linewidth=4.0 if is_analytical else 1.8,
            marker=None if is_analytical else "o",
            markersize=2.8 if not is_analytical else 0.0,
            markevery=max(len(x) // 18, 1) if not is_analytical else None,
            alpha=0.55 if is_analytical else 0.98,
            zorder=draw_order((label, x, y)),
        )
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    out = output_dir / file_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def main() -> int:
    args = parse_common_args("Post-process surface-informed polymer layer verification")
    # This verification intentionally renders each of the four variants when VisIt is
    # available because continuum/CZ and tension/compression are different field tests.
    args.visit_all_variants = True
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    rows_out: list[dict] = []
    summaries: list[dict] = []
    stress_series: list[tuple[str, list[float], list[float]]] = []
    domain_stress_series: list[tuple[str, list[float], list[float]]] = []
    plastic_series: list[tuple[str, list[float], list[float]]] = []
    visit_tex: list[str] = []

    # Put analytical curves first so the comparison is visible even if a run has no CSV yet.
    # Restrict these curves to variants present in the manifest.  This keeps a partial rerun, such
    # as cohesive-only debugging, from plotting unrelated continuum analytical curves and forcing an
    # apparently inconsistent y-axis scale.
    manifest_variant_names = {str(subcase.get("name", "")) for subcase in manifest.get("subcases", [])}
    variants_for_analytical = [variant for variant in VARIANTS if not manifest_variant_names or variant["name"] in manifest_variant_names]
    for variant in variants_for_analytical:
        xs, ss, ps = analytical_series(str(variant["model_kind"]), float(variant["final_global_strain"]))
        stress_series.append(("analytical " + variant["label"], xs, ss))
        plastic_series.append(("analytical " + variant["label"], xs, ps))

    for subcase in manifest.get("subcases", []):
        variant_name = str(subcase.get("name", ""))
        variant = next((v for v in VARIANTS if v["name"] == variant_name), None)
        if variant is None:
            continue
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": variant_name, "label": subcase.get("label"), "error": "run directory not found"})
            continue

        final_strain = float(variant["final_global_strain"])
        reaction_series = read_reaction_stress_y(run_dir, str(variant["model_kind"]), final_strain)
        domain_reaction_series = read_reaction_stress_domain_stretch_y(run_dir, str(variant["model_kind"]), final_strain)
        box_path = find_first(run_dir, ["boxAverageHistory.csv"])
        box_rows = read_csv_numeric(box_path) if box_path else []

        measured_x: list[float] = []
        measured_stress: list[float] = []
        stress_expected: list[float] = []
        stress_errors: list[float] = []
        domain_stretch_x: list[float] = []
        domain_measured_stress: list[float] = []
        plastic_x: list[float] = []
        measured_plastic: list[float] = []
        plastic_expected: list[float] = []
        plastic_errors: list[float] = []

        # Both branches are now loaded through bars/blocks, so boundary reaction
        # stress is the comparable traction diagnostic.  Box averages remain useful
        # for polymer state variables and as a stress fallback for incomplete runs.
        raw_rows: list[dict] = []
        use_reaction_stress = bool(reaction_series)
        stress_source = "reactionHistory" if use_reaction_stress else "boxAverageHistory"
        reaction_source = reaction_measure(run_dir)

        if use_reaction_stress:
            for time, film_strain, stress in reaction_series:
                _film_strain, expected_stress, expected_plastic = expected_at_film_strain(str(variant["model_kind"]), final_strain, film_strain)
                measured_x.append(film_strain)
                measured_stress.append(stress)
                stress_expected.append(expected_stress)
                raw_rows.append(
                    {
                        "variant": variant_name,
                        "label": variant["label"],
                        "time": time,
                        "film_strain_expected": film_strain,
                        "measured_sigma_y": stress,
                        "expected_sigma_y": expected_stress,
                        "stress_error": stress - expected_stress,
                        "measured_equivalent_plastic_strain": "",
                        "expected_equivalent_plastic_strain": expected_plastic,
                        "plastic_strain_error": "",
                        "stress_source": stress_source,
                        "reaction_measure": reaction_source,
                    }
                )
        else:
            for row in box_rows:
                time, film_strain = film_strain_from_history(str(variant["model_kind"]), row, final_strain)
                _film_strain, expected_stress, expected_plastic = expected_at_film_strain(str(variant["model_kind"]), final_strain, film_strain)
                sy = pick_box_material_stress_y(row)
                if sy is None:
                    continue
                measured_x.append(film_strain)
                measured_stress.append(sy)
                stress_expected.append(expected_stress)
                raw_rows.append(
                    {
                        "variant": variant_name,
                        "label": variant["label"],
                        "time": time,
                        "film_strain_expected": film_strain,
                        "measured_sigma_y": sy,
                        "expected_sigma_y": expected_stress,
                        "stress_error": sy - expected_stress,
                        "measured_equivalent_plastic_strain": "",
                        "expected_equivalent_plastic_strain": expected_plastic,
                        "plastic_strain_error": "",
                        "stress_source": stress_source,
                        "reaction_measure": reaction_source,
                    }
                )

        measured_stress, measured_scale = normalize_measured_stress_units(measured_stress, stress_expected)
        stress_errors = [m - e for m, e in zip(measured_stress, stress_expected)]
        if domain_reaction_series:
            domain_stretch_x = [item[1] for item in domain_reaction_series]
            domain_raw_stress = [item[2] for item in domain_reaction_series]
            domain_measured_stress, _domain_scale = normalize_measured_stress_units(domain_raw_stress, [])
        for i, row in enumerate(raw_rows):
            if i < len(measured_stress):
                row["measured_sigma_y"] = measured_stress[i]
                row["stress_error"] = stress_errors[i]
                row["measured_stress_scale_to_mpa"] = measured_scale
                rows_out.append(row)

        status, complete, final_measured_film_strain, target_film_strain = series_completion_status(
            str(variant["model_kind"]),
            final_strain,
            measured_x,
        )

        for row in box_rows:
            time, film_strain = film_strain_from_history(str(variant["model_kind"]), row, final_strain)
            _film_strain, _expected_stress, expected_plastic = expected_at_film_strain(str(variant["model_kind"]), final_strain, film_strain)
            ep = pick_plastic(row)
            if ep is None:
                continue
            plastic_x.append(film_strain)
            measured_plastic.append(ep)
            plastic_expected.append(expected_plastic)
            plastic_errors.append(ep - expected_plastic)

        measured_label_suffix = " measured" if complete else " measured (incomplete)"
        if measured_x:
            stress_series.append((variant["label"] + measured_label_suffix, measured_x, measured_stress))
        if domain_stretch_x:
            domain_stress_series.append((variant["label"] + measured_label_suffix, domain_stretch_x, domain_measured_stress))
        if plastic_x:
            plastic_series.append((variant["label"] + measured_label_suffix, plastic_x, measured_plastic))

        summary = {
            "variant": variant_name,
            "label": variant["label"],
            "status": status,
            "complete": complete,
            "stress_source": stress_source,
            "reaction_measure": reaction_source,
            "target_film_strain": target_film_strain,
            "final_measured_film_strain": final_measured_film_strain,
            "num_stress_samples": len(stress_errors),
            "rms_stress_error": math.sqrt(sum(e * e for e in stress_errors) / len(stress_errors)) if complete and stress_errors else None,
            "max_abs_stress_error": max((abs(e) for e in stress_errors), default=None) if complete else None,
            "final_measured_sigma_y": measured_stress[-1] if measured_stress else None,
            "final_expected_sigma_y": stress_expected[-1] if stress_expected else None,
            "num_plastic_samples": len(plastic_errors),
            "rms_plastic_strain_error": math.sqrt(sum(e * e for e in plastic_errors) / len(plastic_errors)) if complete and plastic_errors else None,
            "max_abs_plastic_strain_error": max((abs(e) for e in plastic_errors), default=None) if complete else None,
            "final_measured_plastic_strain": measured_plastic[-1] if measured_plastic else None,
            "final_expected_plastic_strain": plastic_expected[-1] if plastic_expected else None,
        }
        summaries.append(summary)

        stress_frames = render_visit_frames(
            args,
            run_dir,
            output_dir,
            str(variant["case_name"]),
            "particleStressYY",
            states="initial,quarter,middle,final",
            view="xy",
            colortable="hot_desaturated",
            range_mode="auto",
        )
        visit_tex.extend(
            visit_frames_tex(
                stress_frames,
                output_dir,
                r"Loading-direction stress $\sigma_{yy}$ at initial, quarter, middle, and final states for " + latex_escape(variant["label"]) + ".",
                width="0.23\\linewidth",
                max_frames=4,
            )
        )
        plastic_frames = render_visit_frames(
            args,
            run_dir,
            output_dir,
            str(variant["case_name"]),
            "equivalentPlasticStrain",
            states="initial,quarter,middle,final",
            view="xy",
            colortable="hot_desaturated",
            range_mode="auto",
        )
        visit_tex.extend(
            visit_frames_tex(
                plastic_frames,
                output_dir,
                "Equivalent plastic strain at initial, quarter, middle, and final states for " + latex_escape(variant["label"]) + ".",
                width="0.23\\linewidth",
                max_frames=4,
            )
        )

    plot_metric(output_dir, "surface_polymer_stress_strain.png", "Thin-layer stress-strain response", "film strain", "sigma_y [MPa]", stress_series)
    plot_metric(
        output_dir,
        "surface_polymer_stress_strain_continuum.png",
        "Continuum layer stress-strain response",
        "film strain",
        "sigma_y [MPa]",
        [series for series in stress_series if "continuum" in series[0].lower()],
    )
    plot_metric(
        output_dir,
        "surface_polymer_stress_strain_cohesive.png",
        "Cohesive-zone layer stress-strain response",
        "film strain",
        "sigma_y [MPa]",
        [series for series in stress_series if "cohesive" in series[0].lower()],
    )
    plot_metric(output_dir, "surface_polymer_boundary_stress_domain_stretch.png", "Boundary stress vs domain stretch", "domain stretch F11", "boundary sigma_y [MPa]", domain_stress_series)
    plot_metric(output_dir, "surface_polymer_plastic_strain.png", "Equivalent plastic strain", "film strain", "equivalent plastic strain", plastic_series)
    write_rows(output_dir / "surface_polymer_layer_metrics.csv", rows_out)
    analytical_model = "generic FKM/Viton material-card thin-film update"
    write_json(
        output_dir / "surface_polymer_layer_summary.json",
        {
            "summaries": summaries,
            "film_thickness": FILM_THICKNESS,
            "continuum_gage_length": CONTINUUM_GAGE_LENGTH,
            "cohesive_gage_length": COHESIVE_GAGE_LENGTH,
            "analytical_model": analytical_model,
            "analytical_material": ANALYTICAL_MATERIAL.get("name", "vitonFKM75SurfacePolymer"),
            "analytical_parameters": {
                "bulk_modulus": POLYMER_BULK_MODULUS,
                "shear_modulus": POLYMER_SHEAR_MODULUS,
                "yield_strength": YIELD_STRENGTH,
                "softening_magnitude": SOFTENING_MAGNITUDE,
                "strain_hardening_slope": HARDENING_SLOPE,
                "maximum_stretch": MAXIMUM_STRETCH,
                "pressure_asymmetry": PRESSURE_ASYMMETRY,
            },
        },
    )

    analytical_text = (
        r"The analytical comparison uses the generic FKM/Viton material card constants from "
        r"\texttt{pfw\_materials.py}; the continuum expected curve integrates log film "
        r"strain from the finite-strain update, while the cohesive expected curve uses "
        r"normal jump divided by thickness as in the CZ law."
    )
    tex = [
        r"\paragraph{Quantitative result.} The continuum variants place a finite-thickness continuum polymer layer between elastic loading bars; the cohesive-zone variants place a finite-thickness polymer CZ between elastic blocks.  Measured curves use boundary reactions when present and fall back to material-volume-corrected box averages otherwise.  "
        + analytical_text
        + r"  Cases that do not reach the target film strain are marked incomplete and excluded from error norms.  Stresses are reported in MPa.",
        r"\begin{center}{\scriptsize\resizebox{\linewidth}{!}{\begin{tabular}{llllrrrrr}\toprule Variant & status & stress source & reaction measure & samples & final strain & target strain & RMS stress error & final stress \\\midrule",
    ]
    for s in summaries:
        tex.append(
            latex_escape(s.get("label", s.get("variant", "")))
            + " & "
            + latex_escape(s.get("status", ""))
            + " & "
            + latex_escape(s.get("stress_source", ""))
            + " & "
            + latex_escape(s.get("reaction_measure", ""))
            + " & "
            + compact_or_na(s.get("num_stress_samples", 0), 0)
            + " & "
            + compact_or_na(s.get("final_measured_film_strain"))
            + " & "
            + compact_or_na(s.get("target_film_strain"))
            + " & "
            + compact_or_na(s.get("rms_stress_error"))
            + " & "
            + compact_or_na(s.get("final_measured_sigma_y"))
            + r" \\"
        )
    tex.append(r"\bottomrule\end{tabular}}}\end{center}")
    tex.extend(latex_plot_block(output_dir, "surface_polymer_stress_strain.png", "All variants: stress-strain response."))
    tex.extend(latex_plot_block(output_dir, "surface_polymer_stress_strain_continuum.png", "Continuum variants: stress-strain response."))
    tex.extend(latex_plot_block(output_dir, "surface_polymer_stress_strain_cohesive.png", "Cohesive-zone variants: stress-strain response."))
    tex.extend(latex_plot_block(output_dir, "surface_polymer_boundary_stress_domain_stretch.png", "Measured boundary stress versus imposed domain stretch; this plot uses F11 directly and does not convert the x-axis to film strain."))
    tex.extend(latex_plot_block(output_dir, "surface_polymer_plastic_strain.png", "Equivalent plastic strain histories."))
    tex.extend(visit_tex)
    (output_dir / "surface_polymer_layer_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
