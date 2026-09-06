#!/usr/bin/env python3
"""Shared lightweight post-processing for weak-trace verification folders."""
from __future__ import annotations

import csv
from pathlib import Path
import math
import sys

from mpm_vv_simple_post import (
    compact_float,
    latex_escape,
    load_manifest,
    matplotlib,
    parse_common_args,
    plot_metric,
    read_csv_numeric,
    read_tracer,
    render_visit_frames,
    run_dir_from_subcase,
    tracer_files,
    visit_frames_tex,
    write_json,
    write_rows,
)


WEAK_TRACE_VARIANTS = [
    {"name": "singleField", "label": "single-field baseline"},
    {"name": "traceContactGroups", "label": "weak trace, contact groups"},
    {"name": "falseElasticCZ", "label": "stiff elastic CZ reference"},
]


VARIANT_LINE_STYLES = {
    "singleField": {"color": "#009e73", "linestyle": "-", "linewidth": 4.8, "alpha": 0.50, "zorder": 1.0},
    "traceContactGroups": {"color": "#e69f00", "linestyle": "-", "linewidth": 1.15, "zorder": 2.5},
    "falseElasticCZ": {"color": "#0072b2", "linestyle": "--", "linewidth": 1.55, "zorder": 3.0},
}


def variants_for(case_id: str) -> list[dict]:
    return [
        {
            "name": variant["name"],
            "label": variant["label"],
            "case_name": f"{case_id}_{variant['name']}",
        }
        for variant in WEAK_TRACE_VARIANTS
    ]


def _variant_line_style(diag: dict) -> dict:
    return dict(VARIANT_LINE_STYLES.get(str(diag.get("variant", "")), {"linewidth": 1.2}))


def _deduplicated_legend(axes: list) -> tuple[list, list[str]]:
    handles_out: list = []
    labels_out: list[str] = []
    seen: set[str] = set()
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if label in seen:
                continue
            seen.add(label)
            handles_out.append(handle)
            labels_out.append(label)
    return handles_out, labels_out


def _time_value(row: dict) -> float | None:
    for key in ("time", "Time", "t"):
        try:
            return float(row[key])
        except Exception:
            pass
    return None


def _pick_numeric(row: dict, preferred_tokens: tuple[str, ...]) -> tuple[str, float] | None:
    for key, value in row.items():
        lower = key.lower()
        if all(token in lower for token in preferred_tokens):
            try:
                return key, float(value)
            except Exception:
                pass
    return None


def _final_numeric(row: dict, preferred: tuple[tuple[str, ...], ...], *, allow_fallback: bool = True) -> tuple[str, float] | None:
    for tokens in preferred:
        picked = _pick_numeric(row, tokens)
        if picked is not None:
            return picked
    if not allow_fallback:
        return None
    for key, value in row.items():
        if key.lower() in {"time", "t", "cycle"}:
            continue
        try:
            v = float(value)
        except Exception:
            continue
        if math.isfinite(v):
            return key, v
    return None


def _mean_tracer_series(
    run_dir: Path,
    tracer_prefix: str,
    preferred: tuple[tuple[str, ...], ...],
) -> tuple[str, list[float], list[float], int] | None:
    files = tracer_files(run_dir, tracer_prefix)
    if not files:
        return None

    sums: dict[float, float] = {}
    counts: dict[float, int] = {}
    column = ""
    used_files = 0
    for path in files:
        rows = read_tracer(path)
        file_contributed = False
        for row in rows:
            t = _time_value(row)
            if t is None:
                continue
            picked = _final_numeric(row, preferred, allow_fallback=False)
            if picked is None:
                continue
            key, value = picked
            if not column:
                column = key
            sums[t] = sums.get(t, 0.0) + value
            counts[t] = counts.get(t, 0) + 1
            file_contributed = True
        if file_contributed:
            used_files += 1

    if not counts:
        return None

    times = sorted(counts)
    values = [sums[t] / counts[t] for t in times]
    return f"mean {column} over {used_files} tracers", times, values, used_files


def _dynamic_wave_parameters(config: dict) -> dict:
    left_density = float(config["left_density"])
    right_density = float(config["right_density"])
    left_E = float(config["left_E"])
    right_E = float(config["right_E"])
    poisson = float(config["poisson"])
    left_modulus = left_E * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    right_modulus = right_E * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    left_wave_speed = math.sqrt(left_modulus / left_density)
    right_wave_speed = math.sqrt(right_modulus / right_density)
    left_impedance = left_density * left_wave_speed
    right_impedance = right_density * right_wave_speed
    impedance_sum = left_impedance + right_impedance
    out = dict(config)
    out.update(
        {
            "left_plane_strain_modulus": left_modulus,
            "right_plane_strain_modulus": right_modulus,
            "left_wave_speed": left_wave_speed,
            "right_wave_speed": right_wave_speed,
            "left_impedance": left_impedance,
            "right_impedance": right_impedance,
            "velocity_reflection_coefficient": (left_impedance - right_impedance) / impedance_sum,
            "velocity_transmission_coefficient": 2.0 * left_impedance / impedance_sum,
        }
    )
    return out


def _dynamic_wave_initial_velocity(x: float, params: dict) -> float:
    width = float(params["pulse_width"])
    if width <= 0.0:
        return 0.0
    center = float(params["pulse_center"])
    amplitude = float(params["pulse_velocity"])
    return amplitude * math.exp(-((x - center) / width) ** 2)


def _dynamic_wave_boundary_velocity(time: float, params: dict) -> float:
    if time < 0.0:
        return 0.0
    width = float(params["pulse_time_width"])
    if width <= 0.0:
        return 0.0
    center = float(params["pulse_time_center"])
    amplitude = float(params["pulse_velocity"])
    return amplitude * math.exp(-((time - center) / width) ** 2)


def _dynamic_wave_expected_value(t: float, x: float, component: str, params: dict) -> float:
    interface_x = float(params.get("interface_x", 0.0))
    left_wave_speed = float(params["left_wave_speed"])
    right_wave_speed = float(params["right_wave_speed"])
    left_impedance = float(params["left_impedance"])
    right_impedance = float(params["right_impedance"])
    velocity_reflection = float(params["velocity_reflection_coefficient"])
    velocity_transmission = float(params["velocity_transmission_coefficient"])

    if params.get("source") == "left_boundary_gaussian_velocity":
        source_x = float(params["source_x"])
        if x < interface_x:
            incident_time = t - (x - source_x) / left_wave_speed
            incident = _dynamic_wave_boundary_velocity(incident_time, params)
            reflected = 0.0
            interface_time = t - (interface_x - x) / left_wave_speed
            if interface_time >= 0.0:
                interface_incident_time = interface_time - (interface_x - source_x) / left_wave_speed
                reflected = velocity_reflection * _dynamic_wave_boundary_velocity(interface_incident_time, params)
            if component == "stressXX":
                return -left_impedance * incident + left_impedance * reflected
            return incident + reflected

        interface_time = t - (x - interface_x) / right_wave_speed
        transmitted = 0.0
        if interface_time >= 0.0:
            interface_incident_time = interface_time - (interface_x - source_x) / left_wave_speed
            transmitted = velocity_transmission * _dynamic_wave_boundary_velocity(interface_incident_time, params)
        if component == "stressXX":
            return -right_impedance * transmitted
        return transmitted

    if x < interface_x:
        right_going = 0.5 * _dynamic_wave_initial_velocity(x - left_wave_speed * t, params)
        left_going = 0.5 * _dynamic_wave_initial_velocity(x + left_wave_speed * t, params)
        reflected = 0.0
        reflection_time = t + (x - interface_x) / left_wave_speed
        if reflection_time >= 0.0:
            incident_at_interface = 0.5 * _dynamic_wave_initial_velocity(interface_x - left_wave_speed * reflection_time, params)
            reflected = velocity_reflection * incident_at_interface
        if component == "stressXX":
            return left_impedance * (left_going + reflected - right_going)
        return left_going + right_going + reflected

    transmission_time = t - (x - interface_x) / right_wave_speed
    transmitted = 0.0
    if transmission_time >= 0.0:
        incident_at_interface = 0.5 * _dynamic_wave_initial_velocity(interface_x - left_wave_speed * transmission_time, params)
        transmitted = velocity_transmission * incident_at_interface
    if component == "stressXX":
        return -right_impedance * transmitted
    return transmitted


def _dynamic_wave_tracer_diagnostics(
    run_dir: Path,
    tracer_prefix: str,
    variant: str,
    label: str,
    params: dict,
) -> dict:
    files = tracer_files(run_dir, tracer_prefix)
    out = {"variant": variant, "label": label, "tracers": [], "rows": [], "max_errors": {}}
    for index, path in enumerate(files):
        rows = read_tracer(path)
        if not rows:
            continue
        x0 = float(rows[0]["x"])
        tracer = {"index": index, "name": path.stem, "x0": x0, "series": {}}
        component_series = {
            "velocityX": {"times": [], "values": [], "expected": [], "errors": []},
            "stressXX": {"times": [], "values": [], "expected": [], "errors": []},
        }
        for row in rows:
            t = _time_value(row)
            if t is None:
                continue
            metric_row = {
                "variant": variant,
                "label": label,
                "metric": "dynamic_wave_tracer",
                "tracer": index,
                "tracer_x0": x0,
                "time": t,
            }
            has_value = False
            for component, series in component_series.items():
                if component not in row:
                    continue
                try:
                    measured = float(row[component])
                except Exception:
                    continue
                if not math.isfinite(measured):
                    continue
                expected = _dynamic_wave_expected_value(t, x0, component, params)
                error = measured - expected
                series["times"].append(t)
                series["values"].append(measured)
                series["expected"].append(expected)
                series["errors"].append(error)
                metric_row[component] = measured
                metric_row[f"expected_{component}"] = expected
                metric_row[f"error_{component}"] = error
                has_value = True
            if has_value:
                out["rows"].append(metric_row)
        for component, series in component_series.items():
            if series["times"]:
                tracer["series"][component] = series
                errors = [abs(v) for v in series["errors"]]
                key = f"max_abs_{component}_analytical_error"
                out["max_errors"][key] = max(out["max_errors"].get(key, 0.0), max(errors))
        if tracer["series"]:
            out["tracers"].append(tracer)
    return out


def _dynamic_wave_mean_analytical_series(run_dir: Path, tracer_prefix: str, component: str, params: dict) -> tuple[str, list[float], list[float]] | None:
    files = tracer_files(run_dir, tracer_prefix)
    if not files:
        return None
    sums: dict[float, float] = {}
    counts: dict[float, int] = {}
    used_files = 0
    for path in files:
        rows = read_tracer(path)
        if not rows:
            continue
        x0 = float(rows[0]["x"])
        file_contributed = False
        for row in rows:
            t = _time_value(row)
            if t is None:
                continue
            value = _dynamic_wave_expected_value(t, x0, component, params)
            sums[t] = sums.get(t, 0.0) + value
            counts[t] = counts.get(t, 0) + 1
            file_contributed = True
        if file_contributed:
            used_files += 1
    if not counts:
        return None
    times = sorted(counts)
    values = [sums[t] / counts[t] for t in times]
    return f"analytical mean {component} over {used_files} tracers", times, values


def _plot_dynamic_wave_tracer_component(
    output_dir: Path,
    file_name: str,
    title: str,
    diagnostics: list[dict],
    component: str,
    ylabel: str,
) -> str | None:
    max_tracers = max((len(diag.get("tracers", [])) for diag in diagnostics), default=0)
    if max_tracers <= 0:
        return None

    plt = matplotlib()
    ncols = 2 if max_tracers > 1 else 1
    nrows = int(math.ceil(max_tracers / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.4, 2.7 * nrows), squeeze=False)
    flat_axes = [ax for row in axes for ax in row]
    for tracer_index in range(max_tracers):
        ax = flat_axes[tracer_index]
        plotted_analytical = False
        x0 = None
        for diag in diagnostics:
            tracers_data = diag.get("tracers", [])
            if tracer_index >= len(tracers_data):
                continue
            tracer = tracers_data[tracer_index]
            series = tracer.get("series", {}).get(component)
            if not series:
                continue
            times = series["times"]
            values = series["values"]
            expected = series["expected"]
            x0 = tracer.get("x0", x0)
            ax.plot(times, values, label=str(diag.get("label", diag.get("variant", ""))), **_variant_line_style(diag))
            if not plotted_analytical:
                ax.plot(times, expected, color="k", linestyle="-", linewidth=0.85, label="analytical", zorder=4.0)
                plotted_analytical = True
        ax.set_title(f"tracer {tracer_index}, x0={compact_float(x0, 4)}")
        ax.set_xlabel("time")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)
    for ax in flat_axes[max_tracers:]:
        ax.axis("off")
    handles, labels = _deduplicated_legend(flat_axes)
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)), fontsize=8)
        fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94))
    else:
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.suptitle(title, y=0.995)
    out = output_dir / file_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _read_profile(path: Path) -> dict | None:
    if not path.is_file():
        return None

    rows: list[dict] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            return None
        if len(header) < 2 or header[0].strip().lower() != "time":
            return None
        try:
            x_values = [float(value) for value in header[1:]]
        except Exception:
            return None
        for raw in reader:
            if len(raw) < len(x_values) + 1:
                continue
            try:
                time = float(raw[0])
                values = [float(value) for value in raw[1 : len(x_values) + 1]]
            except Exception:
                continue
            rows.append({"time": time, "values": values})

    if not rows:
        return None
    return {"x": x_values, "rows": rows}


def _nearest_profile_row(rows: list[dict], target_time: float) -> dict | None:
    if not rows:
        return None
    return min(rows, key=lambda row: abs(float(row["time"]) - target_time))


def _dynamic_wave_profile_target_times(params: dict) -> list[float]:
    explicit_times = params.get("profile_times")
    if explicit_times:
        return [float(value) for value in explicit_times]

    stop_time = float(params.get("stop_time", 0.0))
    if params.get("source") == "left_boundary_gaussian_velocity":
        source_x = float(params["source_x"])
        interface_x = float(params["interface_x"])
        left_wave_speed = float(params["left_wave_speed"])
        pulse_center = float(params["pulse_time_center"])
        pulse_width = float(params["pulse_time_width"])
        interface_peak_time = pulse_center + max(0.0, interface_x - source_x) / left_wave_speed
        post_interface_step = max(2.0 * pulse_width, 0.10 * max(interface_peak_time - pulse_center, pulse_width))
        candidates = [
            pulse_center + 0.30 * (interface_peak_time - pulse_center),
            pulse_center + 0.65 * (interface_peak_time - pulse_center),
            interface_peak_time,
            interface_peak_time + post_interface_step,
            interface_peak_time + 2.0 * post_interface_step,
        ]
    else:
        candidates = [0.25 * stop_time, 0.50 * stop_time, 0.75 * stop_time, stop_time]

    out: list[float] = []
    for time in candidates:
        if time < 0.0:
            continue
        if stop_time > 0.0 and time > stop_time:
            continue
        if all(abs(time - previous) > 1.0e-12 for previous in out):
            out.append(time)
    return out


def _dynamic_wave_profile_diagnostics(
    run_dir: Path,
    variant: str,
    label: str,
    direction: str = "x",
) -> dict:
    out = {"variant": variant, "label": label, "profiles": {}}
    for component in ("velocityX", "stressXX"):
        profile = _read_profile(run_dir / f"profile_{direction}_{component}.csv")
        if profile is not None:
            out["profiles"][component] = profile
    return out


def _dynamic_wave_profile_error_summary(profile_diagnostics: dict, params: dict) -> dict:
    out: dict[str, float] = {}
    target_times = _dynamic_wave_profile_target_times(params)
    for component, profile in profile_diagnostics.get("profiles", {}).items():
        max_abs_error = 0.0
        used = False
        for target_time in target_times:
            row = _nearest_profile_row(profile.get("rows", []), target_time)
            if row is None:
                continue
            for x, value in zip(profile.get("x", []), row.get("values", [])):
                if not math.isfinite(x) or not math.isfinite(value):
                    continue
                expected = _dynamic_wave_expected_value(float(row["time"]), x, component, params)
                max_abs_error = max(max_abs_error, abs(value - expected))
                used = True
        if used:
            out[f"profile_{component}_max_abs_selected_time_error"] = max_abs_error
    return out


def _select_profile_times(diagnostics: list[dict], component: str, params: dict) -> list[float]:
    target_times = _dynamic_wave_profile_target_times(params)
    for diag in diagnostics:
        profile = diag.get("profiles", {}).get(component)
        if not profile:
            continue
        rows = profile.get("rows", [])
        selected: list[float] = []
        for target_time in target_times:
            row = _nearest_profile_row(rows, target_time)
            if row is None:
                continue
            actual_time = float(row["time"])
            if all(abs(actual_time - previous) > 1.0e-12 for previous in selected):
                selected.append(actual_time)
        return selected
    return target_times


def _plot_dynamic_wave_profile_component(
    output_dir: Path,
    file_name: str,
    title: str,
    diagnostics: list[dict],
    component: str,
    ylabel: str,
    params: dict,
    error_mode: bool = False,
) -> str | None:
    diagnostics = [diag for diag in diagnostics if diag.get("profiles", {}).get(component)]
    if not diagnostics:
        return None

    selected_times = _select_profile_times(diagnostics, component, params)
    if not selected_times:
        return None

    all_x: list[float] = []
    for diag in diagnostics:
        all_x.extend(diag["profiles"][component]["x"])
    if not all_x:
        return None

    x_min = min(all_x)
    x_max = max(all_x)
    sample_count = 500
    analytical_x = [x_min + (x_max - x_min) * i / float(sample_count - 1) for i in range(sample_count)]

    plt = matplotlib()
    ncols = 2 if len(selected_times) > 1 else 1
    nrows = int(math.ceil(len(selected_times) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.4, 2.9 * nrows), squeeze=False)
    flat_axes = [ax for row in axes for ax in row]
    interface_x = float(params.get("interface_x", 0.0))

    for index, target_time in enumerate(selected_times):
        ax = flat_axes[index]
        for diag in diagnostics:
            profile = diag["profiles"][component]
            row = _nearest_profile_row(profile["rows"], target_time)
            if row is None:
                continue
            x_values = profile["x"]
            y_values = row["values"]
            finite_x: list[float] = []
            finite_y: list[float] = []
            for x, y in zip(x_values, y_values):
                if math.isfinite(x) and math.isfinite(y):
                    finite_x.append(x)
                    if error_mode:
                        finite_y.append(y - _dynamic_wave_expected_value(float(row["time"]), x, component, params))
                    else:
                        finite_y.append(y)
            if finite_x:
                ax.plot(finite_x, finite_y, label=str(diag.get("label", diag.get("variant", ""))), **_variant_line_style(diag))
        if error_mode:
            ax.axhline(0.0, color="k", linestyle="-", linewidth=0.85, label="analytical", zorder=4.0)
        else:
            analytical = [_dynamic_wave_expected_value(target_time, x, component, params) for x in analytical_x]
            ax.plot(analytical_x, analytical, color="k", linestyle="-", linewidth=0.85, label="analytical", zorder=4.0)
        ax.axvline(interface_x, color="0.55", linestyle=":", linewidth=1.0, label="interface", zorder=0.8)
        ax.set_title(f"t = {compact_float(target_time, 5)}")
        ax.set_xlabel("x")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)

    for ax in flat_axes[len(selected_times) :]:
        ax.axis("off")
    handles, labels = _deduplicated_legend(flat_axes)
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)), fontsize=8)
        fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94))
    else:
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.suptitle(title, y=0.995)
    out = output_dir / file_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _lame_from_young_poisson(young: float, poisson: float) -> tuple[float, float, float]:
    shear = young / (2.0 * (1.0 + poisson))
    lame = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    return lame, shear, lame + shear


def _cosine_fraction(time: float, start_time: float, end_time: float) -> float:
    if end_time <= start_time:
        return 1.0
    fraction = min(1.0, max(0.0, (time - start_time) / (end_time - start_time)))
    return 0.5 * (1.0 - math.cos(math.pi * fraction))


def _linear_fraction(time: float, start_time: float, end_time: float) -> float:
    if end_time <= start_time:
        return 1.0
    return min(1.0, max(0.0, (time - start_time) / (end_time - start_time)))


def _ramp_fraction(time: float, params: dict) -> float:
    start_time = float(params.get("load_start_time", 0.0))
    end_time = float(params.get("load_end_time", params.get("stop_time", 1.0)))
    interp = str(params.get("load_interp", "cosine")).strip().lower()
    if interp in {"cosine", "1"}:
        return _cosine_fraction(time, start_time, end_time)
    return _linear_fraction(time, start_time, end_time)


def _radial_pressure_at_time(time: float, params: dict) -> float:
    start_pressure = float(params.get("start_pressure", 0.0))
    end_pressure = float(params.get("external_pressure", params.get("end_pressure", 0.0)))
    return start_pressure + _ramp_fraction(time, params) * (end_pressure - start_pressure)


def _radial_solution_constants(params: dict, pressure: float) -> tuple[float, float, float]:
    a = float(params["inclusion_radius"])
    b = float(params["outer_radius"])
    matrix_E = float(params["matrix_E"])
    inclusion_E = float(params["inclusion_E"])
    poisson = float(params["poisson"])
    _, matrix_shear, matrix_h = _lame_from_young_poisson(matrix_E, poisson)
    _, _, inclusion_h = _lame_from_young_poisson(inclusion_E, poisson)
    denominator = (
        1.0 / (matrix_h * b * b)
        + 1.0 / (matrix_shear * a * a)
        - (1.0 / (b * b) - 1.0 / (a * a)) / inclusion_h
    )
    if denominator == 0.0:
        return -pressure, 0.0, -pressure
    d_matrix = pressure * (1.0 / matrix_h - 1.0 / inclusion_h) / denominator
    c_matrix = -pressure + d_matrix / (b * b)
    c_inclusion = c_matrix - d_matrix / (a * a)
    return c_matrix, d_matrix, c_inclusion


def _radial_expected_components(time: float, x0: float, y0: float, params: dict) -> dict:
    radius = math.hypot(x0, y0)
    a = float(params["inclusion_radius"])
    b = float(params["outer_radius"])
    pressure = _radial_pressure_at_time(time, params)
    c_matrix, d_matrix, c_inclusion = _radial_solution_constants(params, pressure)
    if radius <= a:
        sigma_rr = c_inclusion
        sigma_tt = c_inclusion
    elif radius <= b:
        sigma_rr = c_matrix - d_matrix / (radius * radius)
        sigma_tt = c_matrix + d_matrix / (radius * radius)
    else:
        sigma_rr = math.nan
        sigma_tt = math.nan
    theta = math.atan2(y0, x0) if radius > 0.0 else 0.0
    ct = math.cos(theta)
    st = math.sin(theta)
    return {
        "sigmaRR": sigma_rr,
        "sigmaTT": sigma_tt,
        "stressXX": sigma_rr * ct * ct + sigma_tt * st * st,
        "stressYY": sigma_rr * st * st + sigma_tt * ct * ct,
        "stressXY": (sigma_rr - sigma_tt) * st * ct,
        "pressure": pressure,
    }


def _anti_plane_gamma_at_time(time: float, params: dict) -> float:
    gap = float(params["gap"])
    plate_velocity = float(params["plate_velocity_z"])
    return 2.0 * plate_velocity * max(0.0, time) / gap if gap > 0.0 else 0.0


def _anti_plane_expected_components(time: float, x0: float, y0: float, params: dict) -> dict:
    radius = math.hypot(x0, y0)
    a = float(params["inclusion_radius"])
    matrix_E = float(params["matrix_E"])
    inclusion_E = float(params["inclusion_E"])
    poisson = float(params["poisson"])
    _, matrix_shear, _ = _lame_from_young_poisson(matrix_E, poisson)
    _, inclusion_shear, _ = _lame_from_young_poisson(inclusion_E, poisson)
    gamma = _anti_plane_gamma_at_time(time, params)
    if radius <= a:
        tau_yz = inclusion_shear * gamma * (2.0 * matrix_shear / (matrix_shear + inclusion_shear))
        return {"stressXZ": 0.0, "stressYZ": tau_yz, "gamma": gamma}
    if radius <= 0.0:
        return {"stressXZ": 0.0, "stressYZ": 0.0, "gamma": gamma}
    theta = math.atan2(y0, x0)
    q = (matrix_shear - inclusion_shear) / (matrix_shear + inclusion_shear)
    ar2 = a * a / (radius * radius)
    tau_rz = matrix_shear * gamma * (1.0 - q * ar2) * math.sin(theta)
    tau_tz = matrix_shear * gamma * (1.0 + q * ar2) * math.cos(theta)
    return {
        "stressXZ": tau_rz * math.cos(theta) - tau_tz * math.sin(theta),
        "stressYZ": tau_rz * math.sin(theta) + tau_tz * math.cos(theta),
        "gamma": gamma,
    }


def _stress_components_from_row(row: dict) -> dict | None:
    try:
        return {
            "stressXX": float(row.get("stressXX", 0.0)),
            "stressYY": float(row.get("stressYY", 0.0)),
            "stressZZ": float(row.get("stressZZ", 0.0)),
            "stressYZ": float(row.get("stressYZ", 0.0)),
            "stressXZ": float(row.get("stressXZ", 0.0)),
            "stressXY": float(row.get("stressXY", 0.0)),
        }
    except Exception:
        return None


def _radial_measured_components(row: dict, x0: float, y0: float) -> dict | None:
    stress = _stress_components_from_row(row)
    if stress is None:
        return None
    radius = math.hypot(x0, y0)
    if radius <= 0.0:
        return None
    nx = x0 / radius
    ny = y0 / radius
    tx = -ny
    ty = nx
    sigma_rr = (
        nx * nx * stress["stressXX"]
        + ny * ny * stress["stressYY"]
        + 2.0 * nx * ny * stress["stressXY"]
    )
    sigma_tt = (
        tx * tx * stress["stressXX"]
        + ty * ty * stress["stressYY"]
        + 2.0 * tx * ty * stress["stressXY"]
    )
    out = dict(stress)
    out["sigmaRR"] = sigma_rr
    out["sigmaTT"] = sigma_tt
    return out


def _axis_label(x0: float, y0: float, tolerance: float = 1.0e-10) -> tuple[str, float] | None:
    if abs(y0) <= tolerance and abs(x0) > tolerance:
        return "x", x0
    if abs(x0) <= tolerance and abs(y0) > tolerance:
        return "y", y0
    return None


def _radial_ring_label(radius: float, params: dict, tolerance: float) -> str | None:
    inclusion_radius = float(params["inclusion_radius"])
    outer_radius = float(params["outer_radius"])
    if radius <= inclusion_radius and abs(radius - inclusion_radius) <= tolerance:
        return "interface_inner"
    if radius >= inclusion_radius and abs(radius - inclusion_radius) <= tolerance:
        return "interface_outer"
    if abs(radius - outer_radius) <= tolerance:
        return "outer_surface"
    return None


def _circular_inclusion_tracer_diagnostics(
    run_dir: Path,
    tracer_prefix: str,
    variant: str,
    label: str,
    params: dict,
) -> dict:
    mode = str(params.get("mode", "")).strip().lower()
    axis_tolerance = float(params.get("axis_tolerance", 1.0e-10))
    ring_tolerance = float(params.get("ring_tolerance", 1.5 * axis_tolerance))
    out = {"variant": variant, "label": label, "rows": [], "ring_rows": [], "max_errors": {}}
    quantities = ("sigmaRR", "sigmaTT") if mode == "radial" else ("stressYZ", "stressXZ")
    for index, path in enumerate(tracer_files(run_dir, tracer_prefix)):
        rows = read_tracer(path)
        if not rows:
            continue
        first = rows[0]
        final = rows[-1]
        try:
            x0 = float(first["x"])
            y0 = float(first["y"])
            time = float(final.get("t", final.get("time", 0.0)))
        except Exception:
            continue
        radius = math.hypot(x0, y0)
        axis = _axis_label(x0, y0, axis_tolerance)
        ring = _radial_ring_label(radius, params, ring_tolerance) if mode == "radial" else None
        if axis is None and ring is None:
            continue
        if mode == "radial":
            measured = _radial_measured_components(final, x0, y0)
            expected = _radial_expected_components(time, x0, y0, params)
        else:
            measured = _stress_components_from_row(final)
            expected = _anti_plane_expected_components(time, x0, y0, params)
        if measured is None:
            continue
        row = {
            "variant": variant,
            "label": label,
            "metric": f"circular_{mode}_final_axis_tracer",
            "tracer": path.name,
            "tracerIndex": index,
            "axis": "",
            "coordinate0": "",
            "angleDeg0": math.degrees(math.atan2(y0, x0)) if radius > 0.0 else 0.0,
            "ring": ring or "",
            "radius0": radius,
            "time": time,
            "materialType": final.get("materialType", ""),
        }
        for quantity in quantities:
            value = measured.get(quantity, math.nan)
            target = expected.get(quantity, math.nan)
            row[quantity] = value
            row[f"expected_{quantity}"] = target
            row[f"error_{quantity}"] = value - target if math.isfinite(value) and math.isfinite(target) else math.nan
        if axis is not None:
            axis_name, coordinate = axis
            axis_row = dict(row)
            axis_row["axis"] = axis_name
            axis_row["coordinate0"] = coordinate
            out["rows"].append(axis_row)
        if ring is not None:
            ring_row = dict(row)
            ring_row["metric"] = f"circular_{mode}_final_ring_tracer"
            out["ring_rows"].append(ring_row)
    for quantity in quantities:
        errors = [abs(float(row[f"error_{quantity}"])) for row in out["rows"] if math.isfinite(float(row.get(f"error_{quantity}", math.nan)))]
        if errors:
            out["max_errors"][f"final_{quantity}_max_abs_axis_tracer_error"] = max(errors)
        for ring in ("interface_inner", "interface_outer", "outer_surface"):
            ring_errors = [
                abs(float(row[f"error_{quantity}"]))
                for row in out["ring_rows"]
                if row.get("ring") == ring and math.isfinite(float(row.get(f"error_{quantity}", math.nan)))
            ]
            if ring_errors:
                out["max_errors"][f"final_{quantity}_max_abs_{ring}_tracer_error"] = max(ring_errors)
                out["max_errors"][f"final_{quantity}_rms_{ring}_tracer_error"] = math.sqrt(
                    sum(error * error for error in ring_errors) / float(len(ring_errors))
                )
    return out


def _circular_expected_for_quantity(mode: str, time: float, coordinate: float, axis: str, quantity: str, params: dict) -> float:
    x0 = coordinate if axis == "x" else 0.0
    y0 = coordinate if axis == "y" else 0.0
    if mode == "radial":
        return _radial_expected_components(time, x0, y0, params).get(quantity, math.nan)
    return _anti_plane_expected_components(time, x0, y0, params).get(quantity, math.nan)


def _plot_circular_axis_quantity(
    output_dir: Path,
    file_name: str,
    title: str,
    diagnostics: list[dict],
    params: dict,
    quantity: str,
    ylabel: str,
) -> str | None:
    rows = [row for diag in diagnostics for row in diag.get("rows", []) if quantity in row]
    if not rows:
        return None
    mode = str(params.get("mode", "")).strip().lower()
    axes_to_plot = ["x", "y"]
    plt = matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.7), squeeze=False, sharey=True)
    flat_axes = [ax for row in axes for ax in row]
    outer_radius = float(params.get("outer_radius", max(abs(float(row["coordinate0"])) for row in rows)))
    inclusion_radius = float(params["inclusion_radius"])
    final_time = max(float(row["time"]) for row in rows)
    sample_count = 400

    for ax, axis in zip(flat_axes, axes_to_plot):
        coords = [
            -outer_radius + 2.0 * outer_radius * i / float(sample_count - 1)
            for i in range(sample_count)
        ]
        expected = [_circular_expected_for_quantity(mode, final_time, x, axis, quantity, params) for x in coords]
        ax.plot(coords, expected, color="k", linewidth=0.9, label="analytical", zorder=4.0)
        ax.axvline(-inclusion_radius, color="0.55", linestyle=":", linewidth=1.0, label="interface", zorder=0.8)
        ax.axvline(inclusion_radius, color="0.55", linestyle=":", linewidth=1.0, zorder=0.8)
        for diag in diagnostics:
            pts = [row for row in diag.get("rows", []) if row.get("axis") == axis and math.isfinite(float(row.get(quantity, math.nan)))]
            pts.sort(key=lambda row: float(row["coordinate0"]))
            if not pts:
                continue
            style = _variant_line_style(diag)
            style.update({"marker": "o", "markersize": 3.0})
            ax.plot(
                [float(row["coordinate0"]) for row in pts],
                [float(row[quantity]) for row in pts],
                label=str(diag.get("label", diag.get("variant", ""))),
                **style,
            )
        ax.set_title(f"{axis}-axis, t = {compact_float(final_time, 5)}")
        ax.set_xlabel(f"{axis} coordinate")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)

    handles, labels = _deduplicated_legend(flat_axes)
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)), fontsize=8)
        fig.tight_layout(rect=(0.0, 0.10, 1.0, 0.92))
    else:
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.suptitle(title, y=0.995)
    out = output_dir / file_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _plot_circular_ring_quantity(
    output_dir: Path,
    file_name: str,
    title: str,
    diagnostics: list[dict],
    quantity: str,
    ylabel: str,
) -> str | None:
    ring_rows = [row for diag in diagnostics for row in diag.get("ring_rows", []) if quantity in row]
    if not ring_rows:
        return None
    rings = [
        ("interface_inner", "inside interface"),
        ("interface_outer", "outside interface"),
        ("outer_surface", "outer surface"),
    ]
    rings = [(key, label) for key, label in rings if any(row.get("ring") == key for row in ring_rows)]
    if not rings:
        return None

    plt = matplotlib()
    fig, axes = plt.subplots(len(rings), 1, figsize=(8.0, 2.45 * len(rings)), squeeze=False, sharex=True)
    flat_axes = [ax for row in axes for ax in row]

    for ax, (ring, ring_label) in zip(flat_axes, rings):
        plotted_analytical = False
        for diag in diagnostics:
            pts = [
                row for row in diag.get("ring_rows", [])
                if row.get("ring") == ring and math.isfinite(float(row.get(quantity, math.nan)))
            ]
            pts.sort(key=lambda row: float(row.get("angleDeg0", 0.0)))
            if not pts:
                continue
            angles = [float(row["angleDeg0"]) for row in pts]
            values = [float(row[quantity]) for row in pts]
            expected = [float(row[f"expected_{quantity}"]) for row in pts]
            if not plotted_analytical:
                ax.plot(angles, expected, color="k", linewidth=0.9, label="analytical", zorder=4.0)
                plotted_analytical = True
            style = _variant_line_style(diag)
            style.update({"marker": "o", "markersize": 2.8})
            ax.plot(angles, values, label=str(diag.get("label", diag.get("variant", ""))), **style)
        ax.set_title(ring_label)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)
        ax.set_xlim(-180.0, 180.0)
    flat_axes[-1].set_xlabel("initial polar angle (deg)")

    handles, labels = _deduplicated_legend(flat_axes)
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(4, len(handles)), fontsize=8)
        fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94))
    else:
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.suptitle(title, y=0.995)
    out = output_dir / file_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def postprocess(
    *,
    source_dir: Path,
    default_case_id: str,
    title: str,
    result_prefix: str | None = None,
    description: str | None = None,
    dynamic_wave: dict | None = None,
    circular_inclusion: dict | None = None,
    visit_plots: list[dict] | None = None,
) -> int:
    args = parse_common_args(description or f"Post-process {default_case_id}")
    source_dir = Path(args.source_dir) if args.source_dir else source_dir
    case_id = args.case_id or default_case_id
    output_dir = Path(args.output_dir) if args.output_dir else source_dir / "output" / case_id
    output_dir.mkdir(parents=True, exist_ok=True)
    result_prefix = result_prefix or default_case_id

    variants = variants_for(case_id)
    manifest = load_manifest(source_dir, output_dir, case_id, default_case_id, variants)

    rows_out: list[dict] = []
    summaries: list[dict] = []
    stress_series: list[tuple[str, list[float], list[float]]] = []
    velocity_series: list[tuple[str, list[float], list[float]]] = []
    dynamic_diagnostics: list[dict] = []
    dynamic_profile_diagnostics: list[dict] = []
    dynamic_params = _dynamic_wave_parameters(dynamic_wave) if dynamic_wave else None
    circular_diagnostics: list[dict] = []
    analytical_mean_velocity_added = False
    visit_tex: list[str] = []
    visit_plot_specs = visit_plots if visit_plots is not None else [
        {
            "name": "stress",
            "variable": "Stress",
            "states": "initial,middle,final",
            "view": "xy",
            "colortable": "hot_desaturated",
            "range_mode": "auto",
            "caption": "particle stress smoke-render frames.",
        }
    ]

    for subcase in manifest.get("subcases", []):
        name = str(subcase.get("name", ""))
        label = str(subcase.get("label", name))
        case_name = str(subcase.get("case_name", f"{case_id}_{name}"))
        run_dir = run_dir_from_subcase(source_dir, subcase)
        summary = {
            "variant": name,
            "label": label,
            "run_dir": str(run_dir) if run_dir is not None else "",
            "box_average_samples": 0,
        }
        if run_dir is None:
            summary["status"] = "missing-run-dir"
            summaries.append(summary)
            continue

        box_path = run_dir / "boxAverageHistory.csv"
        data = read_csv_numeric(box_path) if box_path.is_file() else []
        stress_times: list[float] = []
        velocity_times: list[float] = []
        stress_values: list[float] = []
        velocity_values: list[float] = []
        stress_column = ""
        velocity_column = ""
        for row in data:
            t = _time_value(row)
            if t is None:
                continue
            stress = _final_numeric(row, (("stress", "xx"), ("sxx",), ("stress", "xy"), ("sxy",)))
            velocity = _final_numeric(
                row,
                (("velocity", "x"), ("vel", "x"), ("velocity", "y"), ("vel", "y"), ("velocity", "z"), ("vel", "z")),
                allow_fallback=False,
            )
            out_row = {"variant": name, "label": label, "time": t}
            if stress is not None:
                stress_column, stress_value = stress
                out_row["box_average_stress_column"] = stress_column
                out_row["box_average_stress"] = stress_value
                stress_values.append(stress_value)
                stress_times.append(t)
            if velocity is not None:
                velocity_column, velocity_value = velocity
                out_row["velocity_source"] = "boxAverageHistory"
                out_row["velocity_column"] = velocity_column
                out_row["velocity"] = velocity_value
                velocity_values.append(velocity_value)
                velocity_times.append(t)
            if stress is not None or velocity is not None:
                rows_out.append(out_row)

        if not velocity_values:
            tracer_velocity = _mean_tracer_series(
                run_dir,
                f"{result_prefix}Tracer",
                (("velocity", "x"), ("vel", "x"), ("velocity", "y"), ("vel", "y"), ("velocity", "z"), ("vel", "z")),
            )
            if tracer_velocity is not None:
                velocity_column, velocity_times, velocity_values, tracer_count = tracer_velocity
                summary["tracer_velocity_files"] = tracer_count
                for t, value in zip(velocity_times, velocity_values):
                    rows_out.append(
                        {
                            "variant": name,
                            "label": label,
                            "time": t,
                            "velocity_source": "tracer-mean",
                            "velocity_column": velocity_column,
                            "velocity": value,
                        }
                    )
            if dynamic_params is not None and not analytical_mean_velocity_added:
                analytical_mean = _dynamic_wave_mean_analytical_series(run_dir, f"{result_prefix}Tracer", "velocityX", dynamic_params)
                if analytical_mean is not None:
                    velocity_series.append(analytical_mean)
                    analytical_mean_velocity_added = True

        summary["box_average_samples"] = len(data)
        if stress_values:
            summary["stress_column"] = stress_column
            summary["final_box_average_stress"] = stress_values[-1]
            summary["max_abs_box_average_stress"] = max(abs(v) for v in stress_values)
            stress_series.append((label, stress_times, stress_values))
        if velocity_values:
            summary["velocity_column"] = velocity_column
            summary["velocity_source"] = "tracer-mean" if summary.get("tracer_velocity_files") else "boxAverageHistory"
            summary["final_velocity"] = velocity_values[-1]
            summary["max_abs_velocity"] = max(abs(v) for v in velocity_values)
            velocity_series.append((label, velocity_times, velocity_values))
        if dynamic_params is not None:
            diagnostics = _dynamic_wave_tracer_diagnostics(run_dir, f"{result_prefix}Tracer", name, label, dynamic_params)
            if diagnostics.get("tracers"):
                dynamic_diagnostics.append(diagnostics)
                rows_out.extend(diagnostics.get("rows", []))
                summary.update(diagnostics.get("max_errors", {}))
            profile_diagnostics = _dynamic_wave_profile_diagnostics(run_dir, name, label)
            if profile_diagnostics.get("profiles"):
                dynamic_profile_diagnostics.append(profile_diagnostics)
                for component, profile in profile_diagnostics["profiles"].items():
                    summary[f"profile_{component}_samples"] = len(profile.get("rows", []))
                    summary[f"profile_{component}_slices"] = len(profile.get("x", []))
                summary.update(_dynamic_wave_profile_error_summary(profile_diagnostics, dynamic_params))
            else:
                summary["status"] = "missing-profile-output"
                summary["missing_profile_components"] = ["velocityX", "stressXX"]
        if circular_inclusion is not None:
            circular_diag = _circular_inclusion_tracer_diagnostics(
                run_dir,
                f"{result_prefix}Tracer",
                name,
                label,
                circular_inclusion,
            )
            if circular_diag.get("rows"):
                circular_diagnostics.append(circular_diag)
                rows_out.extend(circular_diag.get("rows", []))
                rows_out.extend(circular_diag.get("ring_rows", []))
                summary.update(circular_diag.get("max_errors", {}))

        for spec in visit_plot_specs:
            plot_name = str(spec.get("name", spec.get("variable", "field"))).strip() or "field"
            dedupe_key = spec.get("dedupe_key", f"{name}_{plot_name}")
            frames = render_visit_frames(
                args,
                run_dir,
                output_dir,
                f"{case_name}_{plot_name}",
                str(spec.get("variable", "Stress")),
                states=str(spec.get("states", "initial,middle,final")),
                view=str(spec.get("view", "xy")),
                colortable=str(spec.get("colortable", "hot_desaturated")),
                range_mode=str(spec.get("range_mode", "auto")),
                color_min=spec.get("color_min"),
                color_max=spec.get("color_max"),
                strict_variable=bool(spec.get("strict_variable", False)),
                vector_variables=spec.get("vector_variables"),
                vector_scale=spec.get("vector_scale"),
                vector_count=spec.get("vector_count"),
                vector_stride=spec.get("vector_stride"),
                dedupe_key=None if dedupe_key is None else str(dedupe_key),
            )
            visit_tex.extend(
                visit_frames_tex(
                    frames,
                    output_dir,
                    f"{label} {spec.get('caption', str(spec.get('variable', 'field')) + ' field frames.')}",
                    width=str(spec.get("width", "0.31\\linewidth")),
                    max_frames=int(spec.get("max_frames", 3)),
                )
            )
        summaries.append(summary)

    stress_plot = None
    velocity_plot = None
    dynamic_velocity_plot = None
    dynamic_stress_plot = None
    dynamic_velocity_profile_plot = None
    dynamic_stress_profile_plot = None
    dynamic_velocity_profile_error_plot = None
    dynamic_stress_profile_error_plot = None
    circular_plots: list[str] = []
    circular_ring_plots: list[str] = []
    try:
        stress_plot = plot_metric(output_dir, f"{result_prefix}_box_average_stress.png", title, "time", "box-average stress", stress_series)
        velocity_plot = plot_metric(output_dir, f"{result_prefix}_tracer_velocity.png", title, "time", "mean tracer velocity", velocity_series)
        dynamic_velocity_plot = _plot_dynamic_wave_tracer_component(
            output_dir,
            f"{result_prefix}_tracer_velocity_analytical.png",
            title + ": tracer velocity",
            dynamic_diagnostics,
            "velocityX",
            "velocityX",
        )
        dynamic_stress_plot = _plot_dynamic_wave_tracer_component(
            output_dir,
            f"{result_prefix}_tracer_stress_analytical.png",
            title + ": tracer stress",
            dynamic_diagnostics,
            "stressXX",
            "stressXX",
        )
        if dynamic_params is not None:
            dynamic_velocity_profile_plot = _plot_dynamic_wave_profile_component(
                output_dir,
                f"{result_prefix}_profile_velocity_analytical.png",
                title + ": velocity profiles",
                dynamic_profile_diagnostics,
                "velocityX",
                "velocityX",
                dynamic_params,
            )
            dynamic_stress_profile_plot = _plot_dynamic_wave_profile_component(
                output_dir,
                f"{result_prefix}_profile_stress_analytical.png",
                title + ": stress profiles",
                dynamic_profile_diagnostics,
                "stressXX",
                "stressXX",
                dynamic_params,
            )
            dynamic_velocity_profile_error_plot = _plot_dynamic_wave_profile_component(
                output_dir,
                f"{result_prefix}_profile_velocity_error.png",
                title + ": velocity profile error",
                dynamic_profile_diagnostics,
                "velocityX",
                "velocityX - analytical",
                dynamic_params,
                error_mode=True,
            )
            dynamic_stress_profile_error_plot = _plot_dynamic_wave_profile_component(
                output_dir,
                f"{result_prefix}_profile_stress_error.png",
                title + ": stress profile error",
                dynamic_profile_diagnostics,
                "stressXX",
                "stressXX - analytical",
                dynamic_params,
                error_mode=True,
            )
        if circular_inclusion is not None:
            mode = str(circular_inclusion.get("mode", "")).strip().lower()
            if mode == "radial":
                plot = _plot_circular_axis_quantity(
                    output_dir,
                    f"{result_prefix}_axis_sigma_rr.png",
                    title + ": final radial stress",
                    circular_diagnostics,
                    circular_inclusion,
                    "sigmaRR",
                    "sigma_rr",
                )
                if plot:
                    circular_plots.append(plot)
                plot = _plot_circular_ring_quantity(
                    output_dir,
                    f"{result_prefix}_ring_sigma_rr.png",
                    title + ": final radial stress by ring",
                    circular_diagnostics,
                    "sigmaRR",
                    "sigma_rr",
                )
                if plot:
                    circular_ring_plots.append(plot)
                plot = _plot_circular_axis_quantity(
                    output_dir,
                    f"{result_prefix}_axis_sigma_tt.png",
                    title + ": final hoop stress",
                    circular_diagnostics,
                    circular_inclusion,
                    "sigmaTT",
                    "sigma_tt",
                )
                if plot:
                    circular_plots.append(plot)
                plot = _plot_circular_ring_quantity(
                    output_dir,
                    f"{result_prefix}_ring_sigma_tt.png",
                    title + ": final hoop stress by ring",
                    circular_diagnostics,
                    "sigmaTT",
                    "sigma_tt",
                )
                if plot:
                    circular_ring_plots.append(plot)
            else:
                for quantity in ("stressYZ", "stressXZ"):
                    plot = _plot_circular_axis_quantity(
                        output_dir,
                        f"{result_prefix}_axis_{quantity}.png",
                        title + f": final {quantity}",
                        circular_diagnostics,
                        circular_inclusion,
                        quantity,
                        quantity,
                    )
                    if plot:
                        circular_plots.append(plot)
    except Exception as exc:
        summaries.append({"variant": "postprocess", "status": "plot-failed", "error": str(exc)})

    write_rows(output_dir / f"{result_prefix}_metrics.csv", rows_out)
    summary_doc = {"case_id": case_id, "title": title, "summaries": summaries}
    if dynamic_params is not None:
        analytical_keys = [
                "source",
                "source_x",
                "interface_x",
                "left_density",
                "right_density",
                "left_E",
                "right_E",
                "poisson",
                "pulse_center",
                "pulse_width",
                "pulse_time_center",
                "pulse_time_width",
                "pulse_velocity",
                "stop_time",
                "left_plane_strain_modulus",
                "right_plane_strain_modulus",
                "left_wave_speed",
                "right_wave_speed",
                "left_impedance",
                "right_impedance",
                "velocity_reflection_coefficient",
                "velocity_transmission_coefficient",
        ]
        summary_doc["analytical_solution"] = {key: dynamic_params[key] for key in analytical_keys if key in dynamic_params}
    if circular_inclusion is not None:
        summary_doc["circular_inclusion_reference"] = dict(circular_inclusion)
    write_json(output_dir / f"{result_prefix}_summary.json", summary_doc)

    tex: list[str] = [
        r"\paragraph{Weak-trace diagnostic variants.}",
        "This folder compares the same geometry as a single-field baseline, an explicit-surface weak-trace partition by contact group, and a stiff elastic cohesive-zone reference.",
        r"{\scriptsize\begin{tabular}{llrrr}\toprule Variant & status & box samples & max abs box stress & max abs velocity \\\midrule",
    ]
    for summary in summaries:
        if summary.get("variant") == "postprocess":
            continue
        status = summary.get("status", "available" if summary.get("run_dir") else "missing")
        tex.append(
            latex_escape(summary.get("label", summary.get("variant", "")))
            + " & "
            + latex_escape(status)
            + " & "
            + compact_float(summary.get("box_average_samples", 0), 0)
            + " & "
            + compact_float(summary.get("max_abs_box_average_stress"))
            + " & "
            + compact_float(summary.get("max_abs_velocity"))
            + r" \\"
        )
    tex.append(r"\bottomrule\end{tabular}}")

    if dynamic_params is not None:
        if dynamic_params.get("source") == "left_boundary_gaussian_velocity":
            analytical_text = "Dashed analytical curves use the one-dimensional plane-strain interface solution for the Gaussian velocity pulse imposed at the left boundary by the F-table."
        else:
            analytical_text = "Dashed analytical curves use the one-dimensional plane-strain interface solution for the initialized Gaussian velocity pulse, including the equal initial split into right- and left-going waves."
        tex.extend(
            [
                r"\paragraph{Analytical solution.}",
                analytical_text,
                (
                    r"The velocity coefficients are $R_v="
                    + compact_float(dynamic_params["velocity_reflection_coefficient"])
                    + r"$ and $T_v="
                    + compact_float(dynamic_params["velocity_transmission_coefficient"])
                    + r"$."
                ),
            ]
        )

        profile_error_summaries = [
            summary for summary in summaries
            if "profile_velocityX_max_abs_selected_time_error" in summary or "profile_stressXX_max_abs_selected_time_error" in summary
        ]
        if profile_error_summaries:
            tex.extend(
                [
                    r"\paragraph{Profile residual metric.}",
                    "The table reports the maximum absolute numerical-minus-analytical residual over the selected profile times.",
                    r"{\scriptsize\begin{tabular}{lrr}\toprule Variant & max $|v_x-v_x^{ana}|$ & max $|\sigma_{xx}-\sigma_{xx}^{ana}|$ \\\midrule",
                ]
            )
            for summary in profile_error_summaries:
                tex.append(
                    latex_escape(summary.get("label", summary.get("variant", "")))
                    + " & "
                    + compact_float(summary.get("profile_velocityX_max_abs_selected_time_error"))
                    + " & "
                    + compact_float(summary.get("profile_stressXX_max_abs_selected_time_error"))
                    + r" \\"
                )
            tex.append(r"\bottomrule\end{tabular}}")
    if circular_inclusion is not None:
        mode = str(circular_inclusion.get("mode", "")).strip().lower()
        if mode == "radial":
            tex.extend(
                [
                    r"\paragraph{Circular-inclusion reference.}",
                    "Axis tracer points are compared with the bonded, plane-strain composite-cylinder solution for a circular inclusion of radius "
                    + compact_float(circular_inclusion.get("inclusion_radius"))
                    + " and outer radius "
                    + compact_float(circular_inclusion.get("outer_radius"))
                    + ".  The virtual-fluid pressure event supplies the radial traction used for this first run.",
                    r"{\scriptsize\begin{tabular}{lrr}\toprule Variant & max $|\sigma_{rr}-\sigma_{rr}^{ana}|$ & max $|\sigma_{\theta\theta}-\sigma_{\theta\theta}^{ana}|$ \\\midrule",
                ]
            )
            for summary in summaries:
                if summary.get("variant") == "postprocess":
                    continue
                tex.append(
                    latex_escape(summary.get("label", summary.get("variant", "")))
                    + " & "
                    + compact_float(summary.get("final_sigmaRR_max_abs_axis_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_sigmaTT_max_abs_axis_tracer_error"))
                    + r" \\"
                )
            tex.append(r"\bottomrule\end{tabular}}")
            tex.extend(
                [
                    r"\paragraph{Ring tracer residuals.}",
                    "Ring tracer points are grouped by their initial radius to quantify the local interface and outer-boundary stress error away from the coordinate-axis samples.",
                    r"{\scriptsize\begin{tabular}{lrrrrrr}\toprule Variant & inner $|\sigma_{rr}|_\mathrm{rms}$ & inner $|\sigma_{\theta\theta}|_\mathrm{rms}$ & outer $|\sigma_{rr}|_\mathrm{rms}$ & outer $|\sigma_{\theta\theta}|_\mathrm{rms}$ & boundary $|\sigma_{rr}|_\mathrm{rms}$ & boundary $|\sigma_{\theta\theta}|_\mathrm{rms}$ \\\midrule",
                ]
            )
            for summary in summaries:
                if summary.get("variant") == "postprocess":
                    continue
                tex.append(
                    latex_escape(summary.get("label", summary.get("variant", "")))
                    + " & "
                    + compact_float(summary.get("final_sigmaRR_rms_interface_inner_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_sigmaTT_rms_interface_inner_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_sigmaRR_rms_interface_outer_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_sigmaTT_rms_interface_outer_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_sigmaRR_rms_outer_surface_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_sigmaTT_rms_outer_surface_tracer_error"))
                    + r" \\"
                )
            tex.append(r"\bottomrule\end{tabular}}")
        else:
            tex.extend(
                [
                    r"\paragraph{Circular-inclusion reference.}",
                    "Axis tracer points are compared with the bonded infinite-matrix anti-plane circular-inclusion solution.  This is a local interface diagnostic for the thin 3D slab, not an exact finite-square-boundary solution.",
                    r"{\scriptsize\begin{tabular}{lrr}\toprule Variant & max $|\sigma_{yz}-\sigma_{yz}^{ana}|$ & max $|\sigma_{xz}-\sigma_{xz}^{ana}|$ \\\midrule",
                ]
            )
            for summary in summaries:
                if summary.get("variant") == "postprocess":
                    continue
                tex.append(
                    latex_escape(summary.get("label", summary.get("variant", "")))
                    + " & "
                    + compact_float(summary.get("final_stressYZ_max_abs_axis_tracer_error"))
                    + " & "
                    + compact_float(summary.get("final_stressXZ_max_abs_axis_tracer_error"))
                    + r" \\"
                )
            tex.append(r"\bottomrule\end{tabular}}")

    plots = [p for p in (stress_plot, velocity_plot) if p]
    if plots:
        tex.append(r"\begin{center}")
        for plot_path in plots:
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + Path(plot_path).name + r"}")
        tex.append(r"\end{center}")
    else:
        tex.append(r"\emph{No history plot was produced yet; run the folder variants to populate box-average histories.}")
    dynamic_plots = [p for p in (dynamic_velocity_plot, dynamic_stress_plot) if p]
    if dynamic_plots:
        tex.append(r"\begin{center}")
        for plot_path in dynamic_plots:
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + Path(plot_path).name + r"}")
        tex.append(r"\end{center}")
    profile_plots = [p for p in (dynamic_velocity_profile_plot, dynamic_stress_profile_plot) if p]
    if profile_plots:
        tex.extend(
            [
                r"\paragraph{Wave profiles.}",
                "The profile plots compare Eulerian x-line profiles against the same analytical solution at selected times before and after the incident pulse reaches the material interface.  The single-field baseline is drawn as a thick green underlay, the analytical solution is a thin black foreground curve, and the dotted vertical line marks the interface.",
                r"\begin{center}",
            ]
        )
        for plot_path in profile_plots:
            tex.append(r"\includegraphics[width=0.95\linewidth]{\CaseOutputDir/" + Path(plot_path).name + r"}")
        tex.append(r"\end{center}")
    profile_error_plots = [p for p in (dynamic_velocity_profile_error_plot, dynamic_stress_profile_error_plot) if p]
    if profile_error_plots:
        tex.extend(
            [
                r"\paragraph{Profile residuals.}",
                "The residual plots show numerical minus analytical profile values with the same variant line styles.  The single-field residual is drawn as a thick green underlay and zero residual is a thin black foreground line.",
                r"\begin{center}",
            ]
        )
        for plot_path in profile_error_plots:
            tex.append(r"\includegraphics[width=0.95\linewidth]{\CaseOutputDir/" + Path(plot_path).name + r"}")
        tex.append(r"\end{center}")
    if circular_plots:
        tex.extend(
            [
                r"\paragraph{Axis tracer comparisons.}",
                "These plots show final-time tracer values along the coordinate axes.  The dotted lines mark the two intersections with the circular inclusion interface.",
                r"\begin{center}",
            ]
        )
        for plot_path in circular_plots:
            tex.append(r"\includegraphics[width=0.95\linewidth]{\CaseOutputDir/" + Path(plot_path).name + r"}")
        tex.append(r"\end{center}")
    if circular_ring_plots:
        tex.extend(
            [
                r"\paragraph{Ring tracer comparisons.}",
                "These plots show final-time tracer values around the interface and outer-boundary rings as a function of initial polar angle.",
                r"\begin{center}",
            ]
        )
        for plot_path in circular_ring_plots:
            tex.append(r"\includegraphics[width=0.95\linewidth]{\CaseOutputDir/" + Path(plot_path).name + r"}")
        tex.append(r"\end{center}")
    tex.extend(visit_tex)
    (output_dir / f"{result_prefix}_results.tex").write_text("\n".join(tex) + "\n")
    return 0


__all__ = ["WEAK_TRACE_VARIANTS", "postprocess", "variants_for"]
