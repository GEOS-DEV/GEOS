#!/usr/bin/env python3
"""Aggregate post-processing for the PBC compaction momentum test.

The script reduces each subcase momentum-history CSV to a common table, plots the
total x-momentum versus time, optionally renders VisIt x-velocity contours, and
writes the LaTeX fragment included by test.tex.
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
from pathlib import Path

from variants import *



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Post-process all PBC compaction momentum variants")
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-id", default="pbcCompactionMomentumTest")
    parser.add_argument("--python", dest="python_cmd", default=sys.executable)
    parser.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    parser.add_argument("--no-visit", action="store_true", help="Do not launch VisIt; still include any existing valid velocity PNGs in the LaTeX fragment")
    parser.add_argument("--rerender-visit", action="store_true", help="Re-render VisIt frames even when matching PNGs already exist")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def summarize(history: list[dict[str, object]], force_history: list[dict[str, object]] | None = None) -> dict[str, float]:
    if not history:
        return {"num_rows": 0}
    px = [row_float(row, "particle_momentum_x", 0.0) for row in history]
    times = [row_float(row, "time", float(i)) for i, row in enumerate(history)]
    force_rows_for_summary = force_history or []
    small_fx = [row_float(row, "nodal_internal_force_small_mass_x", 0.0) for row in force_rows_for_summary]
    total_fx = [internal_force_total_x(row) for row in force_rows_for_summary]
    return {
        "num_rows": len(history),
        "initial_time": times[0],
        "final_time": times[-1],
        "initial_particle_momentum_x": px[0],
        "final_particle_momentum_x": px[-1],
        "max_abs_particle_momentum_x": max(abs(v) for v in px),
        "max_abs_small_mass_internal_force_x": max(abs(v) for v in small_fx) if small_fx else 0.0,
        "max_abs_total_internal_force_x": max(abs(v) for v in total_fx) if total_fx else 0.0,
    }


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results = collect_results(source_dir, output_dir, args.case_id)
    write_combined_csv(output_dir, results)
    write_summary(output_dir, results)
    make_plots(output_dir, results)
    run_visit_render(args, source_dir, output_dir, results)
    write_latex(output_dir, results)

    errors = [r for r in results if r.get("error")]
    if errors:
        print("PBC compaction post-processing completed with missing data for: " + ", ".join(r["name"] for r in errors))
    else:
        print("PBC compaction post-processing complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
