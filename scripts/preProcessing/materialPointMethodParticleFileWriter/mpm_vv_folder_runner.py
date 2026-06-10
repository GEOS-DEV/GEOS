#!/usr/bin/env python3
"""Reusable folder-level runner for modern GEOS-MPM verification tests.

A modern verification folder owns one pfw_input file and one runTest script.  The
runTest script passes a list of small variants to :func:`run_folder`; each variant
sets a few environment variables and is dispatched through mpm_vv_case_runner.py.
The folder post-processor is run once after all variants have been prepared or
submitted, so one report section can aggregate all quantitative metrics.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable


def pfw_root(source_dir: Path) -> Path:
    for parent in [source_dir] + list(source_dir.parents):
        if (parent / "particleFileWriter.py").is_file() and (parent / "mpm_vv_case_runner.py").is_file():
            return parent
    raise SystemExit(f"could not find particleFileWriter.py above {source_dir}")


def parse_args(description: str, variants: list[dict]) -> argparse.Namespace:
    choices = [str(v["name"]) for v in variants]
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--case-id", default=None)
    parser.add_argument("--python", dest="python_cmd", default=os.environ.get("PFW_PYTHON", sys.executable))
    parser.add_argument("--force", "-y", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--no-visit", action="store_true")
    parser.add_argument("--post-walltime", default="00:05:00")
    parser.add_argument("--post-partition", default="pdebug")
    parser.add_argument("--walltime", default=None)
    parser.add_argument("--ats-fast", action="store_true")
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("PFW_CASE_JOBS", "1")))
    parser.add_argument("--variant", action="append", choices=choices, help="Run only the selected variant; may be repeated")
    return parser.parse_args()


def variant_case_id(folder_case_id: str, variant: dict) -> str:
    return f"{folder_case_id}__{variant['name']}"


def run_command(cmd: list[str], cwd: Path, env: dict[str, str]) -> tuple[int, str]:
    """Run a child verification command while streaming its output.

    Folder runners dispatch one subcase at a time through mpm_vv_case_runner.py.
    The earlier implementation captured the entire child log with
    subprocess.run(..., stdout=PIPE) and printed it only after the child exited.
    That made a legitimate overwrite prompt from mpm_vv_case_runner invisible to
    the user and made the folder run appear to hang.  Stream the child output one
    character at a time so prompts without trailing newlines are visible, while
    still retaining the full log for the manifest.
    """
    env = dict(env)
    env.setdefault("PYTHONUNBUFFERED", "1")
    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
    )
    output_parts: list[str] = []
    assert proc.stdout is not None
    while True:
        chunk = proc.stdout.read(1)
        if chunk:
            output_parts.append(chunk)
            print(chunk, end="", flush=True)
            continue
        if proc.poll() is not None:
            tail = proc.stdout.read()
            if tail:
                output_parts.append(tail)
                print(tail, end="" if tail.endswith("\n") else "\n", flush=True)
            break
    return proc.wait(), "".join(output_parts)


def read_variant_jobs(source_dir: Path, folder_case_id: str, variant: dict) -> dict:
    cid = variant_case_id(folder_case_id, variant)
    output_dir = source_dir / "output" / cid
    jobs_path = output_dir / f"{variant['case_name']}_jobs.json"
    if not jobs_path.is_file():
        return {}
    try:
        return json.loads(jobs_path.read_text())
    except Exception:
        return {}


def dispatch_variant(args: argparse.Namespace, source_dir: Path, runner: Path, input_file: str, variant: dict) -> dict:
    cid = variant_case_id(args.case_id, variant)
    env = os.environ.copy()
    env.update({str(k): str(v) for k, v in variant.get("env", {}).items()})
    env.setdefault("VV_VARIANT_LABEL", str(variant.get("label", variant["name"])))
    env.setdefault("VV_CASE_NAME", str(variant["case_name"]))

    cmd = [
        args.python_cmd,
        str(runner),
        "--suite", args.suite,
        "--case-id", cid,
        "--input", input_file,
        "--source-dir", str(source_dir),
        "--output-prefix", str(variant["case_name"]),
        "--python", args.python_cmd,
        "--post-walltime", args.post_walltime,
        "--post-partition", args.post_partition,
        "--no-post",
    ]
    if args.force:
        cmd.append("--force")
    if args.prepare_only:
        cmd.append("--prepare-only")
    if args.no_visit:
        cmd.append("--no-visit")
    if args.ats_fast:
        cmd.append("--ats-fast")
    if args.walltime:
        cmd += ["--walltime", args.walltime]

    start_time = time.time()
    print(f"[{args.case_id}] starting {variant['name']} with command: {' '.join(cmd)}", flush=True)
    code, output = run_command(cmd, source_dir, env)
    elapsed_seconds = round(time.time() - start_time, 3)
    print(f"[{args.case_id}] {variant['name']} returncode={code} dispatchElapsed={elapsed_seconds}s", flush=True)
    # The child log has already been streamed by run_command; keep the captured
    # text only for the manifest excerpt.
    jobs = read_variant_jobs(source_dir, args.case_id, variant)
    return {
        "name": variant["name"],
        "label": variant.get("label", variant["name"]),
        "case_id": cid,
        "case_name": variant["case_name"],
        "returncode": code,
        "jobs": jobs,
        "dispatch_elapsed_seconds": elapsed_seconds,
        "log_excerpt": output[-12000:],
    }


def post_command(args: argparse.Namespace, source_dir: Path, output_dir: Path, post_script: str) -> list[str]:
    cmd = [
        args.python_cmd,
        str(source_dir / post_script),
        "--suite", args.suite,
        "--source-dir", str(source_dir),
        "--output-dir", str(output_dir),
        "--case-id", args.case_id,
        "--python", args.python_cmd,
    ]
    if args.no_visit:
        cmd.append("--no-visit")
    return cmd


def submit_or_run_post(args: argparse.Namespace, source_dir: Path, output_dir: Path, manifest: dict, post_script: str, post_job_name: str) -> tuple[str, int]:
    if args.prepare_only:
        return "", 0
    geos_jobs = [str(s.get("jobs", {}).get("geos_job_id", "")) for s in manifest.get("subcases", [])]
    geos_jobs = [j for j in geos_jobs if j]
    cmd = post_command(args, source_dir, output_dir, post_script)
    bank = next((str(s.get("jobs", {}).get("bank", "")) for s in manifest.get("subcases", []) if s.get("jobs", {}).get("bank")), os.environ.get("GEOS_BANK", ""))

    if geos_jobs and shutil.which("sbatch"):
        script = output_dir / f"post_{post_job_name}.sbatch"
        lines = ["#!/bin/bash"]
        if bank:
            lines.append(f"#SBATCH -A {bank}")
        lines += [
            f"#SBATCH -p {args.post_partition}",
            "#SBATCH -N 1",
            "#SBATCH -n 1",
            f"#SBATCH -t {args.post_walltime}",
            f"#SBATCH -J post_{post_job_name[:14]}",
            "#SBATCH --dependency=afterany:" + ":".join(geos_jobs),
            "set -euo pipefail",
            f"cd {source_dir}",
            " ".join(repr(x) for x in cmd),
            "",
        ]
        script.write_text("\n".join(lines))
        proc = subprocess.run(["sbatch", str(script)], cwd=source_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output_dir / f"{post_job_name}_post_submit.log").write_text(proc.stdout)
        print(proc.stdout, end="", flush=True)
        match = re.search(r"Submitted batch job\s+([0-9]+)", proc.stdout)
        return (match.group(1) if match else ""), proc.returncode

    proc = subprocess.run(cmd, cwd=source_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output_dir / f"{post_job_name}_postProcess.log").write_text(proc.stdout)
    print(proc.stdout, end="", flush=True)
    return "", proc.returncode


def run_folder(
    *,
    source_dir: Path,
    input_file: str,
    post_script: str,
    variants: Iterable[dict],
    default_case_id: str,
    description: str | None = None,
    manifest_name: str | None = None,
) -> int:
    variants = list(variants)
    args = parse_args(description or f"Run {default_case_id}", variants)
    if args.case_id is None:
        args.case_id = default_case_id
    root = pfw_root(source_dir)
    runner = root / "mpm_vv_case_runner.py"
    output_dir = source_dir / "output" / args.case_id
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = [v for v in variants if not args.variant or v["name"] in set(args.variant)]
    if not selected:
        raise SystemExit("no variants selected")

    if args.force:
        for child in output_dir.iterdir():
            if child.is_dir() and not child.is_symlink():
                shutil.rmtree(child)
            else:
                child.unlink(missing_ok=True)

    workers = max(1, min(int(args.jobs or 1), len(selected)))
    print(f"[{args.case_id}] dispatching {len(selected)} variant(s) with {workers} worker(s)", flush=True)
    if workers == 1:
        subcases = [dispatch_variant(args, source_dir, runner, input_file, variant) for variant in selected]
    else:
        subcases = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_to_variant = {pool.submit(dispatch_variant, args, source_dir, runner, input_file, variant): variant for variant in selected}
            for future in as_completed(future_to_variant):
                subcases.append(future.result())
        subcases.sort(key=lambda rec: [v["name"] for v in selected].index(rec["name"]))

    manifest = {
        "case_id": args.case_id,
        "suite": args.suite,
        "source_dir": str(source_dir),
        "input": input_file,
        "subcases": subcases,
    }
    output_manifest = output_dir / (manifest_name or f"{default_case_id}_jobs.json")
    output_manifest.write_text(json.dumps(manifest, indent=2, default=str))
    post_job, post_code = submit_or_run_post(args, source_dir, output_dir, manifest, post_script, default_case_id)
    manifest["post_job_id"] = post_job
    manifest["post_returncode"] = post_code
    output_manifest.write_text(json.dumps(manifest, indent=2, default=str))

    if any(int(s.get("returncode", 0)) != 0 for s in subcases):
        return 1
    return int(post_code or 0)
