#!/usr/bin/env python3
"""Run one folder-owned verification input with the common case runner."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run a single folder-scoped MPM verification test")
    p.add_argument("--input", required=True, help="Folder-local pfw_input_*.py file")
    p.add_argument("--output-prefix", required=True, help="Base name for generated output products")
    p.add_argument("--source-dir", default="", help="Folder that owns the verification input; defaults to cwd")
    p.add_argument("--suite", default="verification")
    p.add_argument("--case-id", required=True)
    p.add_argument("--python", dest="python_cmd", default=os.environ.get("PFW_PYTHON", sys.executable))
    p.add_argument("--force", "-y", action="store_true")
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--post-walltime", default="00:05:00")
    p.add_argument("--post-partition", default="pdebug")
    p.add_argument("--walltime", default=None)
    p.add_argument("--ats-fast", action="store_true")
    p.add_argument("--jobs", default="1", help="Accepted for suite compatibility; ignored for single-input folders")
    return p.parse_args()


def find_pfw_root(source_dir: Path) -> Path:
    for parent in [source_dir] + list(source_dir.parents):
        if (parent / "particleFileWriter.py").is_file() and (parent / "mpm_vv_case_runner.py").is_file():
            return parent
    raise SystemExit(f"could not find PFW root above {source_dir}")


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve() if args.source_dir else Path.cwd().resolve()
    runner = find_pfw_root(source_dir) / "mpm_vv_case_runner.py"
    cmd = [
        args.python_cmd,
        str(runner),
        "--suite",
        args.suite,
        "--case-id",
        args.case_id,
        "--input",
        args.input,
        "--source-dir",
        str(source_dir),
        "--output-prefix",
        args.output_prefix,
        "--python",
        args.python_cmd,
        "--post-walltime",
        args.post_walltime,
        "--post-partition",
        args.post_partition,
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
    proc = subprocess.run(cmd, cwd=source_dir)
    return int(proc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
