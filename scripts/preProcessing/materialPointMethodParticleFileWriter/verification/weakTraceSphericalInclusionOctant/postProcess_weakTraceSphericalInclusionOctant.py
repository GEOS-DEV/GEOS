#!/usr/bin/env python3
"""Tracer comparison for the weakTraceSphericalInclusionOctant pressure test."""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path


def _float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_reference(path: Path) -> dict[str, list[dict[str, float | str]]]:
    data: dict[str, list[dict[str, float | str]]] = {}
    for row in _read_csv(path):
        region = row["region"]
        data.setdefault(region, []).append(
            {
                "radius": float(row["radius"]),
                "region": region,
                "displacement_r": float(row["displacement_r"]),
                "stress_rr": float(row["stress_rr"]),
                "stress_tt": float(row["stress_tt"]),
            }
        )
    for rows in data.values():
        rows.sort(key=lambda item: float(item["radius"]))
    return data


def _interpolate_reference(reference: dict[str, list[dict[str, float | str]]], radius: float, material: int) -> dict[str, float | str]:
    region = "inclusion" if material == 1 else "matrix"
    rows = reference[region]
    if radius <= float(rows[0]["radius"]):
        return rows[0]
    if radius >= float(rows[-1]["radius"]):
        return rows[-1]
    for left, right in zip(rows, rows[1:]):
        r0 = float(left["radius"])
        r1 = float(right["radius"])
        if r0 <= radius <= r1:
            if r1 == r0:
                return left
            t = (radius - r0) / (r1 - r0)
            return {
                "radius": radius,
                "region": region,
                "displacement_r": float(left["displacement_r"]) + t * (float(right["displacement_r"]) - float(left["displacement_r"])),
                "stress_rr": float(left["stress_rr"]) + t * (float(right["stress_rr"]) - float(left["stress_rr"])),
                "stress_tt": float(left["stress_tt"]) + t * (float(right["stress_tt"]) - float(left["stress_tt"])),
            }
    raise RuntimeError(f"Could not interpolate reference at r={radius}")


def _stress_components(row: dict[str, str]) -> tuple[float, float, float]:
    x = _float(row, "x")
    y = _float(row, "y")
    z = _float(row, "z")
    radius = math.sqrt(x * x + y * y + z * z)
    if radius > 0.0:
        n = (x / radius, y / radius, z / radius)
    else:
        n = (1.0, 0.0, 0.0)
    sxx = _float(row, "stressXX")
    syy = _float(row, "stressYY")
    szz = _float(row, "stressZZ")
    sxy = _float(row, "stressXY")
    sxz = _float(row, "stressXZ")
    syz = _float(row, "stressYZ")
    sigma_rr = (
        n[0] * (sxx * n[0] + sxy * n[1] + sxz * n[2])
        + n[1] * (sxy * n[0] + syy * n[1] + syz * n[2])
        + n[2] * (sxz * n[0] + syz * n[1] + szz * n[2])
    )
    sigma_tt = (sxx + syy + szz - sigma_rr) / 2.0
    return radius, sigma_rr, sigma_tt


def _parse_run(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        label, path = spec.split("=", 1)
    else:
        path = spec
        label = Path(path).name
    return label, Path(path).expanduser().resolve()


def _runs_from_manifest(output_dir: Path) -> list[tuple[str, Path]]:
    manifests = sorted(output_dir.glob("*_jobs.json"))
    for path in manifests:
        try:
            manifest = json.loads(path.read_text())
        except Exception:
            continue
        runs: list[tuple[str, Path]] = []
        for subcase in manifest.get("subcases", []):
            jobs = subcase.get("jobs", {})
            run_dir = jobs.get("run_dir")
            if run_dir:
                runs.append((str(subcase.get("name", Path(run_dir).name)), Path(run_dir).expanduser().resolve()))
        if runs:
            return runs
    return []


def _default_reference(reference_arg: Path | None, runs: list[tuple[str, Path]], source_dir: Path, output_dir: Path) -> Path:
    if reference_arg is not None:
        return reference_arg.expanduser().resolve()
    if runs:
        candidate = runs[0][1] / "spherical_static_reference.csv"
        if candidate.is_file():
            return candidate
    generated = output_dir / "spherical_static_reference.csv"
    script = source_dir / "spherical_reference.py"
    if script.is_file():
        subprocess.run([sys.executable, str(script), "--output", str(generated)], check=True)
    return generated


def _collect_run(label: str, run_dir: Path, prefix: str, reference: dict[str, list[dict[str, float | str]]]) -> list[dict[str, float | int | str]]:
    rows_out: list[dict[str, float | int | str]] = []
    for path in sorted(run_dir.glob(f"{prefix}_*.csv")):
        rows = _read_csv(path)
        if not rows:
            continue
        initial = rows[0]
        final = rows[-1]
        material = int(round(_float(final, "materialType")))
        radius, sigma_rr, sigma_tt = _stress_components(final)
        expected = _interpolate_reference(reference, radius, material)
        tracer = int(path.stem.split("_")[-1])
        rows_out.append(
            {
                "run": label,
                "tracer": tracer,
                "time": _float(final, "t"),
                "materialType": material,
                "region": str(expected["region"]),
                "initialRadius": math.sqrt(_float(initial, "x") ** 2 + _float(initial, "y") ** 2 + _float(initial, "z") ** 2),
                "finalRadius": radius,
                "x": _float(final, "x"),
                "y": _float(final, "y"),
                "z": _float(final, "z"),
                "sigmaRR": sigma_rr,
                "sigmaTT": sigma_tt,
                "expectedSigmaRR": float(expected["stress_rr"]),
                "expectedSigmaTT": float(expected["stress_tt"]),
                "errorSigmaRR": sigma_rr - float(expected["stress_rr"]),
                "errorSigmaTT": sigma_tt - float(expected["stress_tt"]),
            }
        )
    return rows_out


def _summarize(rows: list[dict[str, float | int | str]]) -> dict[str, object]:
    summary: dict[str, object] = {"runs": {}}
    runs = sorted({str(row["run"]) for row in rows})
    for run in runs:
        run_rows = [row for row in rows if row["run"] == run]
        run_summary: dict[str, object] = {
            "tracers": len(run_rows),
            "time_min": min((float(row["time"]) for row in run_rows), default=math.nan),
            "time_max": max((float(row["time"]) for row in run_rows), default=math.nan),
            "by_material": {},
        }
        for material in sorted({int(row["materialType"]) for row in run_rows}):
            subset = [row for row in run_rows if int(row["materialType"]) == material]
            mat_summary: dict[str, float | int] = {"count": len(subset)}
            for quantity in ("SigmaRR", "SigmaTT"):
                errors = [float(row[f"error{quantity}"]) for row in subset]
                mat_summary[f"max_abs_error_{quantity}"] = max((abs(value) for value in errors), default=math.nan)
                mat_summary[f"rms_error_{quantity}"] = math.sqrt(sum(value * value for value in errors) / len(errors)) if errors else math.nan
            run_summary["by_material"][str(material)] = mat_summary
        summary["runs"][run] = run_summary
    return summary


def _write_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    keys = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _plot(output_dir: Path, rows: list[dict[str, float | int | str]], reference: dict[str, list[dict[str, float | str]]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    styles = {
        "inclusion": {"color": "#0072b2", "marker": "o"},
        "matrix": {"color": "#d55e00", "marker": "s"},
    }
    for quantity, expected_key, ylabel in (
        ("SigmaRR", "stress_rr", "sigma_rr"),
        ("SigmaTT", "stress_tt", "sigma_tt"),
    ):
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharex=True)
        for region, ref_rows in reference.items():
            style = styles.get(region, {})
            xs = [float(row["radius"]) for row in ref_rows]
            ys = [float(row[expected_key]) for row in ref_rows]
            axes[0].plot(xs, ys, color=style.get("color"), linewidth=1.2, label=f"{region} reference")
        for run in sorted({str(row["run"]) for row in rows}):
            for region in ("inclusion", "matrix"):
                subset = [row for row in rows if row["run"] == run and row["region"] == region]
                if not subset:
                    continue
                style = styles.get(region, {})
                axes[0].scatter(
                    [float(row["finalRadius"]) for row in subset],
                    [float(row[f"sigma{quantity[-2:]}"]) for row in subset],
                    marker=style.get("marker", "o"),
                    s=24,
                    label=f"{run} {region}",
                    alpha=0.85,
                )
                axes[1].scatter(
                    [float(row["finalRadius"]) for row in subset],
                    [float(row[f"error{quantity}"]) for row in subset],
                    marker=style.get("marker", "o"),
                    s=24,
                    label=f"{run} {region}",
                    alpha=0.85,
                )
        axes[0].set_ylabel(ylabel)
        axes[1].set_ylabel(f"{ylabel} error")
        for ax in axes:
            ax.set_xlabel("final tracer radius")
            ax.grid(True, alpha=0.25)
        axes[1].axhline(0.0, color="black", linewidth=0.8)
        handles, labels = axes[0].get_legend_handles_labels()
        axes[0].legend(handles, labels, fontsize="x-small")
        fig.tight_layout()
        fig.savefig(output_dir / f"spherical_octant_{ylabel}.png", dpi=180)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", default=[], help="Run as label=/path/to/run or /path/to/run. May be repeated.")
    parser.add_argument("--run-dir", type=Path, help="Single run directory, used if --run is omitted.")
    parser.add_argument("--label", default=None, help="Label for --run-dir.")
    parser.add_argument("--reference", type=Path, help="spherical_static_reference.csv. Defaults to the first run directory.")
    parser.add_argument("--tracer-prefix", default="sphericalOctantTracer")
    parser.add_argument("--output-dir", type=Path, default=Path("spherical_octant_post"))
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--case-id", default="weakTraceSphericalInclusionOctant")
    parser.add_argument("--python", dest="python_cmd", default=sys.executable)
    parser.add_argument("--no-visit", action="store_true", help="Accepted for mpm_vv_folder_runner compatibility.")
    args = parser.parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_specs = list(args.run)
    if not run_specs:
        if args.run_dir is not None:
            run_specs = [f"{args.label or args.run_dir.name}={args.run_dir}"]
        else:
            runs = _runs_from_manifest(output_dir)
            if not runs:
                parser.error("provide --run label=path, --run-dir, or an output directory containing a folder-runner *_jobs.json manifest")
    if run_specs:
        runs = [_parse_run(spec) for spec in run_specs]
    source_dir = args.source_dir.expanduser().resolve()
    reference_path = _default_reference(args.reference, runs, source_dir, output_dir)
    reference = _read_reference(reference_path)

    rows: list[dict[str, float | int | str]] = []
    missing_runs: list[dict[str, str]] = []
    for label, run_dir in runs:
        run_rows = _collect_run(label, run_dir, args.tracer_prefix, reference)
        if not run_rows:
            missing_runs.append({"label": label, "run_dir": str(run_dir), "pattern": f"{args.tracer_prefix}_*.csv"})
        rows.extend(run_rows)
    if missing_runs:
        message = "Missing tracer output for requested run(s): " + ", ".join(f"{item['label']} at {item['run_dir']}" for item in missing_runs)
        if not rows:
            raise SystemExit(message)
        print("WARNING:", message, file=sys.stderr)
    summary = _summarize(rows)
    summary["missing_runs"] = missing_runs
    summary["reference"] = str(reference_path)
    summary["note"] = "Tracer stresses are projected using final tracer position and compared to the static small-strain spherical solution."

    _write_rows(output_dir / "spherical_octant_tracer_stress_errors.csv", rows)
    (output_dir / "spherical_octant_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if rows:
        _plot(output_dir, rows, reference)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 1 if missing_runs else 0


if __name__ == "__main__":
    raise SystemExit(main())
