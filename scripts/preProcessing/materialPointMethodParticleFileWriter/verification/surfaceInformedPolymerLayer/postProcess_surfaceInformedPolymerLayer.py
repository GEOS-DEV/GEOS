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

from dataclasses import dataclass
import math
import os
from pathlib import Path
import sys

import numpy as np

try:
    from scipy.interpolate import PchipInterpolator
except Exception:  # pragma: no cover - scipy may not be available in every post env
    PchipInterpolator = None

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

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
# Fallback area for reaction histories when an older run does not write length_x/length_z.
# Current inputs write both lengths, so this is only a defensive default.
CROSS_SECTION_AREA = DOMAIN_WIDTH * 0.01
STRESS_TO_MPA = 1000.0
ANALYTICAL_TEMPERATURE_C = float(os.environ.get("SURFACE_POLYMER_ANALYTICAL_TEMPERATURE_C", "21.0"))
ANALYTICAL_CRYSTALLINITY_PCT = float(os.environ.get("SURFACE_POLYMER_ANALYTICAL_CRYSTALLINITY_PCT", "7.4"))
ANALYTICAL_DRIVER_PARAMS = {
    "referenceTemperature_C": 26.5,
    "glassTransitionTemperature_C": 26.5,
    "referenceCrystallinity_pct": 7.4,
    "crystallinityTransitionWidth_C": 8.0,
    "scaleTemperaturePoints_C": [-50.0, 0.0, 21.0, 26.5, 40.0, 50.0, 70.0],
    "scaleValues_normalizedToTg": [
        8.787097658379782,
        3.5026233380404235,
        1.6251925860509422,
        1.0,
        0.22216527329736693,
        0.13986649563437228,
        0.06234303240344235,
    ],
    "E_Tg_MPa": 318.0310439541339,
    "yield_Tg_MPa": 7.4892740548914185,
    "shearSoftening_Tg_MPa": 1.175514632547253,
    "hardening_Tg_MPa": 2.048565569441983,
    "shearSofteningShapeParameter1": 0.20104317478877604,
    "shearSofteningShapeParameter2": 3.149845645991683,
    "E_crystallinityCoeff_perPct": 0.14762292766171375,
    "yield_crystallinityCoeff_perPct": 0.19546864238279185,
    "hardeningScaleExponent": 0.851980616403875,
    "pressureAsymmetryAmplitude": 1.0753197990944798,
    "pressureAsymmetryWidth_C": 4.0,
    "bulkModulus_MPa": 2800.0,
    "shearModulus_MPa": 116.07,
    "maximumStretch": 3.0,
    "compressivePressureCapRatio": 0.10,
    "tensilePressureCapRatio": 0.10,
}
POLYMER_BULK_MODULUS = 0.2628333333333331
POLYMER_SHEAR_MODULUS = 0.005291946308724832
POLYMER_CONSTRAINED_MODULUS = POLYMER_BULK_MODULUS + 4.0 * POLYMER_SHEAR_MODULUS / 3.0
YIELD_STRENGTH = 0.0030
SOFTENING_MAGNITUDE = 0.0030
SOFTENING_R1 = 0.30
SOFTENING_R2 = 1.25
HARDENING_SLOPE = 0.0020
PRESSURE_ASYMMETRY = 0.0
MAXIMUM_STRETCH = 2.60



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


@dataclass
class SurfacePolymerDriverState:
    K_MPa: float
    G_MPa: float
    Y_MPa: float
    Ssoft_MPa: float
    H_MPa: float
    eta: float
    r1: float
    r2: float
    lambda_max: float
    compressive_cap_MPa: float
    tensile_cap_MPa: float


class SurfacePolymerCurveDriver:
    """Report analytical driver copied from the standalone polymer curve script."""

    def __init__(self, params: dict[str, object]):
        self.params = dict(params)

    def thermal_scale(self, T_C: float) -> float:
        xp = np.asarray(self.params["scaleTemperaturePoints_C"], dtype=float)
        yp = np.asarray(self.params["scaleValues_normalizedToTg"], dtype=float)
        log_y = np.log(yp)
        if PchipInterpolator is not None:
            return float(np.exp(PchipInterpolator(xp, log_y, extrapolate=True)(T_C)))
        return float(np.exp(np.interp(T_C, xp, log_y, left=log_y[0], right=log_y[-1])))

    def crystallinity_activation(self, T_C: float) -> float:
        Tg = float(self.params["glassTransitionTemperature_C"])
        w = float(self.params.get("crystallinityTransitionWidth_C", 8.0))
        z = max(-60.0, min((T_C - Tg) / w, 60.0))
        return 1.0 / (1.0 + math.exp(-z))

    def crystallinity_factor(self, T_C: float, Xc_pct: float, beta_key: str) -> float:
        beta = float(self.params.get(beta_key, 0.0))
        Xref = float(self.params["referenceCrystallinity_pct"])
        return max(0.05, 1.0 + beta * (Xc_pct - Xref) * self.crystallinity_activation(T_C))

    def eta_T(self, T_C: float) -> float:
        Tg = float(self.params["glassTransitionTemperature_C"])
        w = float(self.params.get("pressureAsymmetryWidth_C", 4.0))
        eta0 = float(self.params.get("pressureAsymmetryAmplitude", 0.0))
        return eta0 * math.exp(-0.5 * ((T_C - Tg) / w) ** 2)

    def material_state(self, T_C: float, Xc_pct: float) -> SurfacePolymerDriverState:
        S_T = self.thermal_scale(T_C)
        CY = self.crystallinity_factor(T_C, Xc_pct, "yield_crystallinityCoeff_perPct")
        K = float(self.params.get("bulkModulus_MPa", 2800.0))
        G = float(self.params.get("shearModulus_MPa", self.params["E_Tg_MPa"] / (2.0 * (1.0 + 0.37))))
        comp_ratio = float(self.params.get("compressivePressureCapRatio", -1.0))
        tens_ratio = float(self.params.get("tensilePressureCapRatio", -1.0))
        return SurfacePolymerDriverState(
            K_MPa=K,
            G_MPa=G,
            Y_MPa=float(self.params["yield_Tg_MPa"] * S_T * CY),
            Ssoft_MPa=float(self.params["shearSoftening_Tg_MPa"] * S_T),
            H_MPa=float(self.params["hardening_Tg_MPa"] * S_T ** float(self.params["hardeningScaleExponent"])),
            eta=self.eta_T(T_C),
            r1=float(self.params["shearSofteningShapeParameter1"]),
            r2=float(self.params["shearSofteningShapeParameter2"]),
            lambda_max=float(self.params.get("maximumStretch", 3.0)),
            compressive_cap_MPa=(comp_ratio * K if comp_ratio >= 0.0 else float("inf")),
            tensile_cap_MPa=(tens_ratio * K if tens_ratio >= 0.0 else float("inf")),
        )

    @staticmethod
    def chain_stretch_confined(eps: float) -> float:
        J = max(math.exp(eps), 1.0e-300)
        return max(J ** (-1.0 / 3.0), J ** (2.0 / 3.0))

    @staticmethod
    def chain_stretch_cz(eps_n: float, gamma: float = 0.0) -> float:
        F = np.array([[1.0, gamma, 0.0], [0.0, math.exp(eps_n), 0.0], [0.0, 0.0, 1.0]], dtype=float)
        J = max(float(np.linalg.det(F)), 1.0e-300)
        Cbar = (J ** (-2.0 / 3.0)) * (F.T @ F)
        vals = np.linalg.eigvalsh(Cbar)
        return float(math.sqrt(max(vals.max(), 0.0)))

    @staticmethod
    def hardening_measure(lambda_chain: float) -> float:
        return max(0.0, lambda_chain * lambda_chain - 1.0 / max(lambda_chain, 1.0e-30))

    def base_flow_strength(self, kappa: float, lambda_chain: float, material: SurfacePolymerDriverState) -> float:
        if material.r1 <= 0.0:
            soft = 0.0
        else:
            soft = material.Ssoft_MPa * math.exp(-((max(kappa, 0.0) / material.r1) ** material.r2))
        return material.Y_MPa + soft + material.H_MPa * self.hardening_measure(lambda_chain)

    def confined_series(
        self,
        model: str,
        final_true_strain: float,
        npts: int,
        T_C: float,
        Xc_pct: float,
    ) -> tuple[list[float], list[float], list[float]]:
        material = self.material_state(T_C, Xc_pct)
        true_strains = np.linspace(0.0, final_true_strain, npts)
        sign = 1.0 if final_true_strain >= 0.0 else -1.0
        stresses: list[float] = []
        kappas: list[float] = []
        engineering_strains: list[float] = []
        kappa_previous = 0.0

        for eps in true_strains:
            lambda_chain = self.chain_stretch_cz(float(eps)) if model == "cz" else self.chain_stretch_confined(float(eps))
            if lambda_chain >= material.lambda_max and stresses:
                stress = float("nan")
                kappa = kappa_previous
            else:
                p_t = material.K_MPa * float(eps)
                p_eff = max(p_t, -material.compressive_cap_MPa) if p_t < 0.0 else min(p_t, material.tensile_cap_MPa)
                eps_eq_total = 2.0 * abs(float(eps)) / 3.0
                base = self.base_flow_strength(kappa_previous, lambda_chain, material)
                q_limit = max(0.0, base - material.eta * p_eff)
                q_trial = 3.0 * material.G_MPa * eps_eq_total
                q = min(q_trial, q_limit)
                kappa = max(kappa_previous, eps_eq_total - q / max(3.0 * material.G_MPa, 1.0e-30))
                stress = p_t + sign * (2.0 / 3.0) * q
                kappa_previous = kappa

            engineering_strains.append(math.exp(float(eps)) - 1.0)
            stresses.append(stress)
            kappas.append(kappa)

        return engineering_strains, stresses, kappas


def softening(kappa: float) -> float:
    ratio = max(kappa, 0.0) / max(SOFTENING_R1, 1.0e-30)
    return SOFTENING_MAGNITUDE * math.exp(-ratio ** SOFTENING_R2)


def stretch_hardening(film_strain: float) -> float:
    lam = max(1.0 + film_strain, 1.0)
    lam = min(lam, MAXIMUM_STRETCH)
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
                     kappa_old: float) -> tuple[list[float], list[float], float]:
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
    flow = YIELD_STRENGTH + softening(kappa_old) + stretch_hardening(total_film_strain)
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
        flow_candidate = YIELD_STRENGTH + softening(kappa_candidate) + stretch_hardening(total_film_strain)
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
        flow = YIELD_STRENGTH + softening(kappa) + stretch_hardening(film_strain)
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

def analytical_series(model_kind: str, final_strain: float, n: int = 251) -> tuple[list[float], list[float], list[float]]:
    final_engineering_strain = imposed_film_strain(model_kind, final_strain)
    final_stretch = max(1.0 + final_engineering_strain, 1.0e-12)
    final_true_strain = math.log(final_stretch)
    driver = SurfacePolymerCurveDriver(ANALYTICAL_DRIVER_PARAMS)
    return driver.confined_series(
        "cz" if model_kind == "cohesive" else "continuum",
        final_true_strain,
        n,
        ANALYTICAL_TEMPERATURE_C,
        ANALYTICAL_CRYSTALLINITY_PCT,
    )


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
    return None


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


def read_reaction_stress_y(run_dir: Path, model_kind: str, final_strain: float) -> list[tuple[float, float, float]]:
    path = find_first(run_dir, ["reactionHistory.csv"])
    if path is None:
        return []
    rows = read_csv_numeric(path)
    out: list[tuple[float, float, float]] = []
    for row in rows:
        clean = {str(k).strip(): v for k, v in row.items()}
        time, film_strain = film_strain_from_history(model_kind, clean, final_strain)
        ry_minus = clean.get("Ry-", clean.get("reactionYMinus", None))
        ry_plus = clean.get("Ry+", clean.get("reactionYPlus", None))
        values = []
        for value in [ry_minus, ry_plus]:
            try:
                values.append(float(value))
            except Exception:
                pass
        if values:
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
            out.append((time, film_strain, signed_force / max(area, 1.0e-30)))
    out.sort(key=lambda item: item[0])
    filtered: list[tuple[float, float, float]] = []
    last_t = -math.inf
    for time, film_strain, stress in out:
        if time > last_t:
            filtered.append((time, film_strain, stress))
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

    Continuum box averages are written in material stress units, while boundary
    reaction histories in this verification workflow are already engineering
    MPa after force/area reduction.  The ratio check keeps the post-processor
    robust when a run provides only one of those two histories.
    """
    pairs = [(abs(m), abs(e)) for m, e in zip(measured, expected_mpa) if abs(m) > 1.0e-14 and abs(e) > 1.0e-8]
    if not pairs:
        return measured, 1.0
    ratios = sorted(m / e for m, e in pairs)
    ratio = ratios[len(ratios) // 2]
    if ratio < 0.02:
        scale = STRESS_TO_MPA
    elif ratio > 50.0:
        scale = 1.0 / STRESS_TO_MPA
    else:
        scale = 1.0
    return [scale * value for value in measured], scale


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
        box_path = find_first(run_dir, ["boxAverageHistory.csv"])
        box_rows = read_csv_numeric(box_path) if box_path else []

        measured_x: list[float] = []
        measured_stress: list[float] = []
        stress_expected: list[float] = []
        stress_errors: list[float] = []
        plastic_x: list[float] = []
        measured_plastic: list[float] = []
        plastic_expected: list[float] = []
        plastic_errors: list[float] = []

        # Both branches are now loaded through bars/blocks, so boundary reaction
        # stress is the comparable traction diagnostic.  Box averages remain useful
        # for polymer state variables and as a stress fallback for incomplete runs.
        raw_rows: list[dict] = []
        use_reaction_stress = bool(reaction_series)

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
                        "stress_source": "reactionHistory",
                    }
                )
        else:
            for row in box_rows:
                time, film_strain = film_strain_from_history(str(variant["model_kind"]), row, final_strain)
                _film_strain, expected_stress, expected_plastic = expected_at_film_strain(str(variant["model_kind"]), final_strain, film_strain)
                sy = pick_stress_y(row)
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
                        "stress_source": "boxAverageHistory",
                    }
                )

        measured_stress, measured_scale = normalize_measured_stress_units(measured_stress, stress_expected)
        stress_errors = [m - e for m, e in zip(measured_stress, stress_expected)]
        for i, row in enumerate(raw_rows):
            if i < len(measured_stress):
                row["measured_sigma_y"] = measured_stress[i]
                row["stress_error"] = stress_errors[i]
                row["measured_stress_scale_to_mpa"] = measured_scale
                rows_out.append(row)

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

        if measured_x:
            stress_series.append((variant["label"] + " measured", measured_x, measured_stress))
        if plastic_x:
            plastic_series.append((variant["label"] + " measured", plastic_x, measured_plastic))

        summary = {
            "variant": variant_name,
            "label": variant["label"],
            "num_stress_samples": len(stress_errors),
            "rms_stress_error": math.sqrt(sum(e * e for e in stress_errors) / len(stress_errors)) if stress_errors else None,
            "max_abs_stress_error": max((abs(e) for e in stress_errors), default=None),
            "final_measured_sigma_y": measured_stress[-1] if measured_stress else None,
            "final_expected_sigma_y": stress_expected[-1] if stress_expected else None,
            "num_plastic_samples": len(plastic_errors),
            "rms_plastic_strain_error": math.sqrt(sum(e * e for e in plastic_errors) / len(plastic_errors)) if plastic_errors else None,
            "max_abs_plastic_strain_error": max((abs(e) for e in plastic_errors), default=None),
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
    plot_metric(output_dir, "surface_polymer_plastic_strain.png", "Equivalent plastic strain", "film strain", "equivalent plastic strain", plastic_series)
    write_rows(output_dir / "surface_polymer_layer_metrics.csv", rows_out)
    write_json(
        output_dir / "surface_polymer_layer_summary.json",
        {
            "summaries": summaries,
            "film_thickness": FILM_THICKNESS,
            "continuum_gage_length": CONTINUUM_GAGE_LENGTH,
            "cohesive_gage_length": COHESIVE_GAGE_LENGTH,
            "analytical_model": "attached surface polymer curve driver, confined path",
            "analytical_temperature_C": ANALYTICAL_TEMPERATURE_C,
            "analytical_crystallinity_pct": ANALYTICAL_CRYSTALLINITY_PCT,
            "analytical_parameters": ANALYTICAL_DRIVER_PARAMS,
        },
    )

    tex = [
        r"\paragraph{Quantitative result.} The continuum variants place a finite-thickness continuum polymer layer between elastic loading bars; the cohesive-zone variants place a finite-thickness polymer CZ between elastic blocks.  Both measured curves use boundary reactions.  The analytical comparison uses the attached surface-polymer curve driver in confined thin-film mode at "
        + compact_float(ANALYTICAL_TEMPERATURE_C)
        + r"$^\circ$C and "
        + compact_float(ANALYTICAL_CRYSTALLINITY_PCT)
        + r"\% crystallinity; the driver integrates true strain and the plots report engineering film strain for comparison with the GEOS deformation history.  Stresses are reported in MPa.",
        r"{\scriptsize\begin{tabular}{lrrrr}\toprule Variant & samples & RMS stress error & max stress error & final stress \\\midrule",
    ]
    for s in summaries:
        tex.append(
            latex_escape(s.get("label", s.get("variant", "")))
            + " & "
            + compact_float(s.get("num_stress_samples", 0), 0)
            + " & "
            + compact_float(s.get("rms_stress_error"))
            + " & "
            + compact_float(s.get("max_abs_stress_error"))
            + " & "
            + compact_float(s.get("final_measured_sigma_y"))
            + r" \\"
        )
    tex.append(r"\bottomrule\end{tabular}}")
    for name in [
        "surface_polymer_stress_strain.png",
        "surface_polymer_stress_strain_continuum.png",
        "surface_polymer_stress_strain_cohesive.png",
        "surface_polymer_plastic_strain.png",
    ]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    tex.extend(visit_tex)
    (output_dir / "surface_polymer_layer_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
