#!/usr/bin/env python3
"""Post-process one GEOS-MPM verification case."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
OPTIONAL_LEGACY_PLOT_WARNING = 'optional legacy plotting script failed'
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Post-process one MPM verification case")
    p.add_argument("--suite", default="verification")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--source-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--job-id", default="")
    p.add_argument("--case-id", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--initial-variable", default="auto")
    p.add_argument("--final-variable", default="auto")
    p.add_argument("--view", default="xy")
    p.add_argument("--mesh", default="CellRegion1")
    p.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    p.add_argument("--python", dest="python_cmd", default=os.environ.get("PFW_PYTHON", "/usr/tce/bin/python3"))
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--skip-render", action="store_true")
    p.add_argument("--copy-only", action="store_true")
    p.add_argument("--skip-slurm-check", action="store_true")
    p.add_argument("--force", "-y", action="store_true")
    return p.parse_args()


def log(args, msg):
    print(f"[{args.case_id} post] {msg}", flush=True)


def warn(args, msg):
    print(f"[{args.case_id} post] WARNING: {msg}", file=sys.stderr, flush=True)


def plot_database_exists(run_dir: Path) -> bool:
    silo = run_dir / "siloFiles"
    return silo.is_dir() and (any(silo.glob("mpm_cpdi_*")) or any(silo.glob("mpm_*")))


def slurm_file_for(run_dir: Path, job_id: str):
    if not job_id:
        return None
    for p in [run_dir / f"slurm-{job_id}.out", run_dir / f"slurm-{job_id}_0.out"]:
        if p.is_file():
            return p
    found = sorted(run_dir.glob(f"slurm-{job_id}*.out"))
    return found[0] if found else None


def check_geos_completed(args, run_dir: Path) -> bool:
    if args.skip_slurm_check or not args.job_id:
        return True
    out = run_dir / f"{args.output_prefix}_geos_slurm_check.log"
    slurm = slurm_file_for(run_dir, args.job_id)
    records = []
    if slurm and slurm.is_file():
        data = slurm.read_text(errors="replace")
        tail = data[-50000:]
        records.append(f"slurm={slurm}")
        fail_pat = r"(XML parsing error|EXCEPTION|ERROR|CANCELLED|TIME LIMIT|OUT_OF_MEMORY|NODE_FAIL|Segmentation fault|Floating point error|MPI_Abort|Killed)"
        success_pat = r"(Job complete|GEOSX executed successfully|problem run complete|Total elapsed time|run complete)"
        failure = re.search(fail_pat, tail, re.I)
        success = re.search(success_pat, tail, re.I)
        records.append("success_marker=" + str(bool(success)))
        records.append("failure_marker=" + (failure.group(0) if failure else ""))
        if failure:
            important = []
            for line in tail.splitlines():
                if re.search(fail_pat + r"|missing required attribute|contains unused attribute|invalid within|Error cause", line, re.I):
                    important.append(line)
            records.append("important_tail:\n" + "\n".join(important[-60:]))
        out.write_text("\n".join(records) + "\n")
        if failure and not success:
            warn(args, "failure marker found in GEOS slurm output; copying available products only")
            return False
        if success:
            return True
        if plot_database_exists(run_dir):
            warn(args, "no success marker found, but Silo plot database exists; continuing")
            return True
    else:
        warn(args, f"slurm output for job {args.job_id} was not found")
        if plot_database_exists(run_dir):
            return True
    if shutil.which("sacct") and args.job_id:
        proc = subprocess.run(["sacct", "-j", args.job_id, "-X", "-n", "-P", "-o", "JobID,State,ExitCode"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        line = next((l for l in proc.stdout.splitlines() if l.strip()), "")
        if line:
            parts = line.split("|")
            state = parts[1] if len(parts) > 1 else ""
            exit_code = parts[2] if len(parts) > 2 else ""
            with out.open("a") as f:
                f.write(f"sacct_state={state}\nsacct_exit={exit_code}\n")
            if state == "COMPLETED" and exit_code == "0:0":
                return True
    return plot_database_exists(run_dir)


def detect_visit(requested: str):
    candidates = [requested, os.environ.get("VISIT_CMD", ""), os.environ.get("VISIT_COMMAND", ""), "/usr/gapps/visit/bin/visit", "visit"]
    for c in candidates:
        if not c:
            continue
        if os.path.isabs(c) and os.access(c, os.X_OK):
            return c
        found = shutil.which(c)
        if found:
            return found
    return None


def choose_variables(args):
    cid = args.case_id.lower()
    src = str(args.source_dir).lower()
    family = cid.split("__", 1)[0]
    initial = args.initial_variable
    final = args.final_variable
    if initial != "auto" and final != "auto":
        return initial, final
    if "geomechanics" in cid or "temperaturetable" in cid:
        ini, fin = "Pressure", "PlasticStrainMagnitude"
    elif "ceramic" in cid or "sizeeffect" in cid:
        ini, fin = "StrengthScale", "Damage"
    elif "cohesive" in cid or "polymercz" in cid:
        ini, fin = "MaterialType", "Damage"
    elif "vonmises" in cid:
        ini, fin = "MaterialType", "PlasticStrainMagnitude"
    elif "gas" in cid or "fluid" in cid:
        ini, fin = "Density", "Pressure"
    elif "temperature" in cid or "polymerheal" in cid or "shrinkage" in cid:
        ini, fin = "Temperature", "Temperature"
    elif "momentum" in cid or "periodic" in cid or "partition" in cid or "contact" in cid:
        ini, fin = "Velocity", "Velocity"
    else:
        ini, fin = "MaterialType", "Damage"
    return (ini if initial == "auto" else initial), (fin if final == "auto" else final)


def run_visit_one(args, run_dir: Path, frame_dir: Path, state: str, variable: str, range_mode: str, logf):
    if args.no_visit or args.skip_render or args.copy_only:
        return
    visit = detect_visit(args.visit_cmd)
    script = run_dir / "pfw_visit_render.py"
    if not visit or not script.is_file():
        warn(args, "VisIt or pfw_visit_render.py not available; skipping render")
        return
    cmd = [
        visit, "-nowin", "-cli", "-s", str(script),
        "--run-dir", str(run_dir),
        "--output-dir", str(frame_dir),
        "--case-name", args.output_prefix,
        "--states", state,
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
        warn(args, f"VisIt exited with {proc.returncode}; continuing to copy any PNGs")


def run_visit(args, run_dir: Path, frame_dir: Path, geos_ok: bool):
    if not geos_ok and not plot_database_exists(run_dir):
        return
    frame_dir.mkdir(parents=True, exist_ok=True)
    initial, final = choose_variables(args)
    log_path = run_dir / f"{args.output_prefix}_visit_render.log"
    with log_path.open("w") as logf:
        run_visit_one(args, run_dir, frame_dir, "initial", initial, "auto", logf)
        final_range = "unit" if "damage" in final.lower() else "auto"
        run_visit_one(args, run_dir, frame_dir, "final", final, final_range, logf)


def run_legacy_plot_scripts(args, run_dir: Path):
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    logs = []
    for script in sorted(run_dir.glob("plot*.py")):
        log_file = run_dir / f"{args.output_prefix}_{script.stem}.log"
        proc = subprocess.run([args.python_cmd, str(script.name)], cwd=run_dir, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_file.write_text(proc.stdout + f"\nreturncode={proc.returncode}\n")
        logs.append(str(log_file))
        if proc.returncode != 0:
            warn(args, f"legacy plot script {script.name} exited with {proc.returncode}")
    return logs


def generate_plots(args, run_dir: Path, source_dir: Path):
    sys.path.insert(0, str(run_dir))
    products = []
    try:
        import mpm_vv_plot_tools as tools
        plan = tools.make_plot_plan(source_dir, args.case_id, args.case_id.split("__", 1)[0])
        tools.write_plan(run_dir / f"{args.output_prefix}_plot_requirements.json", plan)
        products.extend(tools.generate_standard_plots(run_dir, args.output_prefix))
    except Exception as exc:
        warn(args, f"generic plotting failed: {exc}")
    products.extend(run_legacy_plot_scripts(args, run_dir))
    return products


def copy_products(args, run_dir: Path, output_dir: Path, frame_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    patterns = ["*.png", "*.csv", "*.json", "*.log", "*.txt"]
    seen = set()
    for base in [frame_dir, run_dir]:
        if not base.exists():
            continue
        for pattern in patterns:
            for src in sorted(base.glob(pattern)):
                if src.is_file() and src.resolve() not in seen:
                    seen.add(src.resolve())
                    dst = output_dir / src.name
                    shutil.copy2(src, dst)
                    copied.append(dst.name)
    (output_dir / f"{args.output_prefix}_copied_products.json").write_text(json.dumps(copied, indent=2))
    log(args, f"copied {len(copied)} product(s) to {output_dir}")


def main():
    args = parse_args()
    run_dir = Path(os.path.expanduser(os.path.expandvars(args.run_dir))).resolve()
    source_dir = Path(os.path.expanduser(os.path.expandvars(args.source_dir))).resolve()
    output_dir = Path(os.path.expanduser(os.path.expandvars(args.output_dir))).resolve()
    frame_dir = run_dir / "visit_frames"
    geos_ok = check_geos_completed(args, run_dir)
    if not args.copy_only:
        generate_plots(args, run_dir, source_dir)
        run_visit(args, run_dir, frame_dir, geos_ok)
    copy_products(args, run_dir, output_dir, frame_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
