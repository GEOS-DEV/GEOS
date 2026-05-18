#!/usr/bin/env python3
"""Generic robust post-processing for staged MPM examples."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Post-process one staged MPM example")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--job-id", default="")
    p.add_argument("--case-name", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--initial-variable", default="Damage")
    p.add_argument("--final-variable", default="Damage")
    p.add_argument("--view", default="xy")
    p.add_argument("--mesh", default="CellRegion1")
    p.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--skip-render", action="store_true")
    p.add_argument("--copy-only", action="store_true")
    p.add_argument("--skip-slurm-check", action="store_true")
    p.add_argument("--force", "-y", action="store_true")
    return p.parse_args()


def log(args, msg):
    print(f"[{args.case_name} post] {msg}", flush=True)


def warn(args, msg):
    print(f"[{args.case_name} post] WARNING: {msg}", file=sys.stderr, flush=True)


def expand(path):
    return Path(os.path.expandvars(os.path.expanduser(str(path)))).resolve()


def plot_database_exists(run_dir: Path) -> bool:
    silo = run_dir / "siloFiles"
    return silo.is_dir() and (any(silo.glob("mpm_cpdi_*")) or any(silo.glob("mpm_*")))


def slurm_file_for(run_dir: Path, job_id: str) -> Path | None:
    if not job_id:
        return None
    candidates = [run_dir / f"slurm-{job_id}.out", run_dir / f"slurm-{job_id}_0.out"]
    for c in candidates:
        if c.is_file():
            return c
    found = sorted(run_dir.glob(f"slurm-{job_id}*.out"))
    return found[0] if found else None


def check_geos_completed(args, run_dir: Path) -> bool:
    if args.skip_slurm_check or not args.job_id:
        return True
    log_path = run_dir / f"{args.output_prefix}_geos_slurm_check.log"
    slurm = slurm_file_for(run_dir, args.job_id)
    records = []
    if slurm and slurm.is_file():
        log(args, "checking GEOS slurm output: " + str(slurm))
        data = slurm.read_text(errors="replace")
        tail = data[-30000:]
        records.append("slurm=" + str(slurm))
        failure = re.search(r"(CANCELLED|TIME LIMIT|OUT_OF_MEMORY|NODE_FAIL|Segmentation fault|Traceback|EXCEPTION|MPI_Abort)", tail, re.I)
        success = re.search(r"(Job complete|GEOSX executed successfully|problem run complete|Total elapsed time|run complete)", tail, re.I)
        records.append("success_marker=" + str(bool(success)))
        records.append("failure_marker=" + (failure.group(0) if failure else ""))
        log_path.write_text("\n".join(records) + "\n")
        if failure and not success:
            warn(args, "failure marker found in GEOS slurm output; post-processing will copy existing products only")
            return False
        if success:
            return True
        if plot_database_exists(run_dir):
            warn(args, "no success marker found, but Silo plot database exists; continuing")
            return True
    else:
        warn(args, f"slurm output for job {args.job_id} was not found")
        if plot_database_exists(run_dir):
            warn(args, "Silo plot database exists; continuing")
            return True
    if shutil.which("sacct") and args.job_id:
        proc = subprocess.run(["sacct", "-j", args.job_id, "-X", "-n", "-P", "-o", "JobID,State,ExitCode"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        line = next((l for l in proc.stdout.splitlines() if l.strip()), "")
        if line:
            parts = line.split("|")
            state = parts[1] if len(parts) > 1 else ""
            exit_code = parts[2] if len(parts) > 2 else ""
            with log_path.open("a") as f:
                f.write(f"sacct_state={state}\nsacct_exit={exit_code}\n")
            log(args, f"sacct state={state} exit={exit_code}")
            if state == "COMPLETED" and exit_code == "0:0":
                return True
            if state in {"RUNNING", "PENDING", "CONFIGURING", "COMPLETING"}:
                return False
    return plot_database_exists(run_dir)


def detect_visit(requested: str) -> str | None:
    candidates = []
    if requested:
        candidates.append(requested)
    candidates += [os.environ.get("VISIT_CMD", ""), os.environ.get("VISIT_COMMAND", ""), "/usr/gapps/visit/bin/visit", "visit"]
    for c in candidates:
        if not c:
            continue
        if os.path.isabs(c) and os.access(c, os.X_OK):
            return c
        found = shutil.which(c)
        if found:
            return found
    return None


def run_visit_command(args, run_dir: Path, frame_dir: Path, states: str, variable: str, range_mode: str, logf) -> None:
    custom = run_dir / f"visitRender_{args.output_prefix}.py"
    if custom.is_file():
        visit = detect_visit(args.visit_cmd)
        if not visit:
            warn(args, "VisIt command not found; skipping render")
            return
        frame_dir.mkdir(parents=True, exist_ok=True)
        log_file = run_dir / f"{args.output_prefix}_visit_render.log"
        cmd = [visit, "-nowin", "-cli", "-s", str(custom), "--run-dir", str(run_dir), "--output-dir", str(frame_dir), "--case-name", args.output_prefix, "--states", "initial,final", "--mesh", args.mesh, "--list-databases"]
        log(args, "running VisIt: " + " ".join(cmd))
        with log_file.open("w") as logf:
            logf.write("Running: " + " ".join(cmd) + "\n")
            proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            logf.write(proc.stdout)
            print(proc.stdout, end="")
            if proc.returncode != 0:
                warn(args, f"VisIt command exited with {proc.returncode}; continuing to copy any generated PNGs")
        return
    script = run_dir / "pfw_visit_render.py"
    visit = detect_visit(args.visit_cmd)
    if not script.is_file():
        warn(args, "pfw_visit_render.py not found in run directory; skipping VisIt")
        return
    if not visit:
        warn(args, "VisIt command not found; skipping render")
        return
    cmd = [
        visit, "-nowin", "-cli", "-s", str(script),
        "--run-dir", str(run_dir),
        "--output-dir", str(frame_dir),
        "--case-name", args.output_prefix,
        "--states", states,
        "--view", args.view,
        "--variable", variable,
        "--mesh", args.mesh,
        "--colortable", "hot_desaturated",
        "--range-mode", range_mode,
        "--list-databases",
    ]
    log(args, "running VisIt: " + " ".join(cmd))
    logf.write("Running: " + " ".join(cmd) + "\n")
    proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    logf.write(proc.stdout)
    print(proc.stdout, end="")
    if proc.returncode != 0:
        warn(args, f"VisIt command exited with {proc.returncode}; continuing to copy any generated PNGs")


def run_visit(args, run_dir: Path, frame_dir: Path) -> None:
    if args.no_visit or args.skip_render or args.copy_only:
        return
    frame_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / f"{args.output_prefix}_visit_render.log"
    with log_file.open("w") as logf:
        run_visit_command(args, run_dir, frame_dir, "initial", args.initial_variable, "auto", logf)
        final_range = "unit" if "damage" in args.final_variable.lower() else "auto"
        run_visit_command(args, run_dir, frame_dir, "final", args.final_variable, final_range, logf)


def parse_float(value: str):
    try:
        return float(value)
    except Exception:
        return None


def plot_reactions(args, run_dir: Path) -> Path | None:
    path = run_dir / "reactionHistory.csv"
    if not path.is_file():
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        warn(args, f"matplotlib unavailable for reaction plot: {exc}")
        return None
    with path.open(newline="") as f:
        reader = csv.reader(f)
        rows = [r for r in reader if r]
    if len(rows) < 2:
        return None
    header = [h.strip() for h in rows[0]]
    numeric = []
    for row in rows[1:]:
        vals = [parse_float(x) for x in row]
        if any(v is not None for v in vals):
            numeric.append(vals)
    if not numeric:
        return None
    time_idx = 0
    for i, h in enumerate(header):
        if h.lower() in {"time", "t"} or "time" in h.lower():
            time_idx = i
            break
    candidates = []
    for i, h in enumerate(header):
        low = h.lower()
        if i == time_idx:
            continue
        if "reaction" in low or "force" in low or "ry" in low or "_y" in low or low.endswith("y"):
            candidates.append(i)
    if not candidates:
        candidates = [i for i in range(len(header)) if i != time_idx][:6]
    x = [row[time_idx] if time_idx < len(row) else None for row in numeric]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    plotted = 0
    for i in candidates[:8]:
        y = [row[i] if i < len(row) else None for row in numeric]
        if all(v is None for v in y):
            continue
        ax.plot(x, y, label=header[i] if i < len(header) else f"col{i}")
        plotted += 1
    if plotted == 0:
        plt.close(fig)
        return None
    ax.set_xlabel(header[time_idx] if time_idx < len(header) else "time")
    ax.set_ylabel("reaction / force")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = run_dir / f"{args.output_prefix}_reactions.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def copy_products(args, run_dir: Path, output_dir: Path, frame_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    seen_sources = set()
    patterns = [f"{args.output_prefix}*.png", f"{args.output_prefix}*.csv", f"{args.output_prefix}*.log", f"{args.output_prefix}*.json", "*.png", "*.log"]
    for base in [frame_dir, run_dir]:
        if not base.exists():
            continue
        for pattern in patterns:
            for src in sorted(base.glob(pattern)):
                if src.is_file() and src.resolve() not in seen_sources:
                    seen_sources.add(src.resolve())
                    dst = output_dir / src.name
                    shutil.copy2(src, dst)
                    copied.append(dst.name)
    for name in ["reactionHistory.csv", "boxAverageHistory.csv"]:
        src = run_dir / name
        if src.is_file():
            shutil.copy2(src, output_dir / name)
            copied.append(name)
    if copied:
        for name in sorted(set(copied)):
            log(args, "copied " + name)
    else:
        warn(args, "no PNG/CSV/log products were found to copy")
        if frame_dir.exists():
            warn(args, "visit_frames contents: " + ", ".join(p.name for p in frame_dir.iterdir()))


def main():
    args = parse_args()
    run_dir = expand(args.run_dir)
    output_dir = expand(args.output_dir)
    frame_dir = run_dir / "visit_frames"
    log(args, "run directory: " + str(run_dir))
    log(args, "output directory: " + str(output_dir))
    log(args, "visit frames directory: " + str(frame_dir))
    ok = check_geos_completed(args, run_dir)
    if not ok:
        warn(args, "GEOS completion check failed or job is not complete; copying existing products only")
        args.copy_only = True
    if not args.copy_only:
        plot_reactions(args, run_dir)
        run_visit(args, run_dir, frame_dir)
    copy_products(args, run_dir, output_dir, frame_dir)
    payload = {"case": args.case_name, "generated_at": datetime.now().isoformat(timespec="seconds"), "job_id": args.job_id, "run_dir": str(run_dir), "output_dir": str(output_dir)}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{args.output_prefix}_postprocess.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
