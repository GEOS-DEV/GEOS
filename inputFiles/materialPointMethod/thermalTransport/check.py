#!/usr/bin/env python3
"""Check thermal diagnostics and tracer histories after running a generated job."""
import argparse
import csv
import json
import math
from pathlib import Path
import re


def check(folder, log_name="run.log"):
    folder = Path(folder)
    metadata = json.loads((folder / "case.json").read_text())
    case = metadata["case"]
    log = (folder / log_name).read_text()
    diagnostics = re.findall(r"MPMThermal[^\n]*balanceResidual=([-+\deE.]+)", log)
    if not diagnostics:
        raise AssertionError("No MPMThermal energy diagnostics found")
    if any(not math.isfinite(float(value)) or abs(float(value)) > 1e-8 for value in diagnostics):
        raise AssertionError("Thermal energy balance failed")
    paths = sorted(folder.glob("thermalTracer_*.csv"))
    if not paths:
        raise AssertionError("No temperature tracers found")
    first, last = [], []
    for path in paths:
        with path.open() as stream:
            records = list(csv.DictReader(stream))
        if not records:
            raise AssertionError(f"Empty tracer {path}")
        row = {key: float(value) for key, value in records[-1].items()}
        if row["t"] < 0.99 * metadata["end_time"]:
            raise AssertionError("Solver did not reach the requested end time")
        t = row["temperature"]
        if not math.isfinite(t) or t < 0:
            raise AssertionError("Invalid absolute temperature")
        if case == "source" and abs(t - (350 + 4 * row["t"])) > 1e-7:
            raise AssertionError("Particle source does not match analytical heating")
        if case == "insulated" and abs(t - (400 if row["x"] < 0.5 else 200)) > 1e-7:
            raise AssertionError("Insulated interface transferred heat")
        if case == "isothermal" and row["x"] < 0.5 and abs(t - 400) > 1e-8:
            raise AssertionError("Isothermal particle changed temperature")
        if case == "bulk":
            exact = 350 + 50 * math.cos(math.pi * row["x"]) * math.exp(-metadata["diffusivity"] * math.pi ** 2 * row["t"])
            if abs(t - exact) > 3.0:
                raise AssertionError("Bulk diffusion exceeded coarse-grid analytical tolerance")
        first.append(records[0])
        last.append(row)
    if case in ("groups", "dfg", "materials", "isothermal"):
        cold = [row["temperature"] for row in last if row["x"] > 0.5]
        if not cold or max(cold) <= 200 + 1e-6:
            raise AssertionError("Closed stationary contact transferred no heat")
        if any(row["temperature"] < 200 - 1e-5 or row["temperature"] > 400 + 1e-5 for row in last):
            raise AssertionError("Contact temperature overshoot")
    print(f"PASS {case}: {len(diagnostics)} energy balances, {len(paths)} tracers")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job", type=Path)
    parser.add_argument("--log", default="run.log")
    args = parser.parse_args()
    check(args.job, args.log)
