#!/usr/bin/env python3
"""Post-process the surface-informed polymer layer verification.

The folder runs four variants: continuum/CZ and tension/compression.  The
quantitative comparison uses the same one-dimensional film model used by the new
constitutive kernels: a constrained elastic predictor, scalar plastic limiting,
softening as a function of accumulated plastic strain, and stretch hardening.  The
VisIt render requests show the loading-direction stress and equivalent plastic
strain at four output states for each variant when VisIt is available.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

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
STOP_TIME = 1.0
GAGE_LENGTH = 1.0
FILM_THICKNESS = 0.10
CROSS_SECTION_AREA = 0.50 * (1.0 / 32.0)  # width times nominal plane-strain thickness.
POLYMER_BULK_MODULUS = 0.2628333333333331
POLYMER_SHEAR_MODULUS = 0.005291946308724832
POLYMER_CONSTRAINED_MODULUS = POLYMER_BULK_MODULUS + 4.0 * POLYMER_SHEAR_MODULUS / 3.0
YIELD_STRENGTH = 0.0030
SOFTENING_MAGNITUDE = 0.0030
SOFTENING_R1 = 0.30
SOFTENING_R2 = 1.25
HARDENING_SLOPE = 0.0020
MAXIMUM_STRETCH = 2.60


def final_global_strain(variant_name: str) -> float:
    variant = next((v for v in VARIANTS if v["name"] == variant_name), None)
    return float(variant.get("final_global_strain", 0.0)) if variant else 0.0


def imposed_film_strain(global_strain: float) -> float:
    """Reconstruct the film strain from imposed gage strain.

    The elastic blocks are intentionally stiff, so the analytical comparison treats
    the polymer film as receiving the imposed displacement jump.  This makes the
    continuum and cohesive-zone variants use the same local one-dimensional
    reference solution.
    """
    return global_strain * GAGE_LENGTH / FILM_THICKNESS


def softening(kappa: float) -> float:
    ratio = max(kappa, 0.0) / max(SOFTENING_R1, 1.0e-30)
    return SOFTENING_MAGNITUDE * math.exp(-ratio ** SOFTENING_R2)


def stretch_hardening(film_strain: float) -> float:
    lam = max(1.0 + film_strain, 1.0)
    lam = min(lam, MAXIMUM_STRETCH)
    return HARDENING_SLOPE * (lam * lam - 1.0 / lam)


def local_update(film_strain: float, plastic_strain_old: float, kappa_old: float) -> tuple[float, float, float]:
    """One-dimensional constrained update matching the verification material inputs."""
    trial = POLYMER_CONSTRAINED_MODULUS * (film_strain - plastic_strain_old)
    flow = YIELD_STRENGTH + softening(kappa_old) + stretch_hardening(film_strain)
    if abs(trial) <= flow:
        return trial, plastic_strain_old, kappa_old
    sign = 1.0 if trial >= 0.0 else -1.0
    stress = sign * flow
    plastic_strain_new = film_strain - stress / POLYMER_CONSTRAINED_MODULUS
    kappa_new = max(kappa_old, abs(plastic_strain_new))
    # Re-evaluate once because softening depends on the updated plastic strain.
    flow = YIELD_STRENGTH + softening(kappa_new) + stretch_hardening(film_strain)
    stress = sign * min(abs(trial), flow)
    plastic_strain_new = film_strain - stress / POLYMER_CONSTRAINED_MODULUS
    kappa_new = max(kappa_old, abs(plastic_strain_new))
    return stress, plastic_strain_new, kappa_new


def analytical_series(final_strain: float, n: int = 251) -> tuple[list[float], list[float], list[float]]:
    xs: list[float] = []
    stresses: list[float] = []
    kappas: list[float] = []
    plastic_strain = 0.0
    kappa = 0.0
    for i in range(n):
        tfrac = i / max(n - 1, 1)
        film_strain = imposed_film_strain(final_strain * tfrac)
        stress, plastic_strain, kappa = local_update(film_strain, plastic_strain, kappa)
        xs.append(film_strain)
        stresses.append(stress)
        kappas.append(kappa)
    return xs, stresses, kappas


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


def read_reaction_stress_y(run_dir: Path) -> list[tuple[float, float]]:
    path = find_first(run_dir, ["reactionHistory.csv"])
    if path is None:
        return []
    rows = read_csv_numeric(path)
    out: list[tuple[float, float]] = []
    for row in rows:
        clean = {str(k).strip(): v for k, v in row.items()}
        try:
            time = float(clean.get("time", clean.get("Time", 0.0)))
        except Exception:
            continue
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
                area = abs(float(clean.get("length_x", 0.50)) * float(clean.get("length_z", 1.0 / 32.0)))
            except Exception:
                area = CROSS_SECTION_AREA
            out.append((time, signed_force / max(area, CROSS_SECTION_AREA, 1.0e-30)))
    out.sort(key=lambda item: item[0])
    filtered: list[tuple[float, float]] = []
    last_t = -math.inf
    for time, stress in out:
        if time > last_t:
            filtered.append((time, stress))
            last_t = time
    return filtered


def expected_at_time(final_strain: float, time: float) -> tuple[float, float, float]:
    n = 251
    target = max(0.0, min(time / STOP_TIME, 1.0))
    xs, ss, ps = analytical_series(final_strain, n=n)
    idx = min(max(int(round(target * (n - 1))), 0), n - 1)
    return xs[idx], ss[idx], ps[idx]


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
    for variant in VARIANTS:
        xs, ss, ps = analytical_series(float(variant["final_global_strain"]))
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
        reaction_series = read_reaction_stress_y(run_dir)
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

        if reaction_series:
            for time, stress in reaction_series:
                film_strain, expected_stress, expected_plastic = expected_at_time(final_strain, time)
                measured_x.append(film_strain)
                measured_stress.append(stress)
                stress_expected.append(expected_stress)
                stress_errors.append(stress - expected_stress)
                rows_out.append(
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
                    }
                )
        else:
            for row in box_rows:
                time = float(row.get("time", row.get("t", 0.0)))
                film_strain, expected_stress, expected_plastic = expected_at_time(final_strain, time)
                sy = pick_stress_y(row)
                if sy is None:
                    continue
                measured_x.append(film_strain)
                measured_stress.append(sy)
                stress_expected.append(expected_stress)
                stress_errors.append(sy - expected_stress)
                rows_out.append(
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
                    }
                )

        for row in box_rows:
            time = float(row.get("time", row.get("t", 0.0)))
            film_strain, _expected_stress, expected_plastic = expected_at_time(final_strain, time)
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

    plot_metric(output_dir, "surface_polymer_stress_strain.png", "Thin-layer stress-strain response", "film strain", "sigma_y", stress_series)
    plot_metric(output_dir, "surface_polymer_plastic_strain.png", "Equivalent plastic strain", "film strain", "equivalent plastic strain", plastic_series)
    write_rows(output_dir / "surface_polymer_layer_metrics.csv", rows_out)
    write_json(
        output_dir / "surface_polymer_layer_summary.json",
        {
            "summaries": summaries,
            "film_thickness": FILM_THICKNESS,
            "constrained_modulus": POLYMER_CONSTRAINED_MODULUS,
            "yield_strength": YIELD_STRENGTH,
            "softening_magnitude": SOFTENING_MAGNITUDE,
            "strain_hardening_slope": HARDENING_SLOPE,
        },
    )

    tex = [
        r"\paragraph{Quantitative result.} The analytical comparison is a one-dimensional thin-film update with constrained modulus $M=K+4G/3$, scalar plastic limiting, softening, and stretch hardening matching the material inputs.",
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
    for name in ["surface_polymer_stress_strain.png", "surface_polymer_plastic_strain.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    tex.extend(visit_tex)
    (output_dir / "surface_polymer_layer_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
