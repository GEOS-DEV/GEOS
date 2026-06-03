#!/usr/bin/env python3
"""Stage, generate, submit, and post-process one verification case."""
from __future__ import annotations

import argparse
import getpass
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


def log(case, msg):
    print(f"[{case}] {msg}", flush=True)


def fail(case, msg):
    print(f"[{case}] ERROR: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(1)


def expand_path(p):
    return Path(os.path.expandvars(os.path.expanduser(str(p)))).resolve()


def parse_args():
    p = argparse.ArgumentParser(description="Run one MPM verification case")
    p.add_argument("--suite", default="verification")
    p.add_argument("--case-id", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--source-dir", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--python", dest="python_cmd", default=os.environ.get("PFW_PYTHON", "/usr/tce/bin/python3"))
    p.add_argument("--force", "-y", action="store_true")
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--post-walltime", default='00:05:00')
    p.add_argument("--post-partition", default="pdebug")
    p.add_argument("--walltime", default=None, help="Optional GEOS walltime override; omitted means use the pfw_input value")
    p.add_argument("--no-submit", action="store_true")
    p.add_argument("--no-post", action="store_true", help="Do not submit the generic per-case post-processing job")
    p.add_argument("--ats-fast", action="store_true", help="Append explicit ultra-fast ATS sizing overrides to the staged input")
    return p.parse_args()


def find_pfw_root(source_dir: Path) -> Path:
    for p in [source_dir] + list(source_dir.parents):
        if (p / "particleFileWriter.py").is_file() and (p / "pfw_geometryObjects.py").is_file():
            return p
    fail("verification", f"could not find PFW root above {source_dir}")


def load_userdefs(pfw_root: Path):
    candidates = []
    for name in (os.environ.get("USER", ""), os.environ.get("LOGNAME", ""), getpass.getuser()):
        if name:
            candidates.append(pfw_root / f"userDefs_{name}.py")
    candidates += sorted(pfw_root.glob("userDefs_*.py"))
    for path in candidates:
        if path.is_file():
            spec = importlib.util.spec_from_file_location("_pfw_userdefs", path)
            mod = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(mod)
            return path, mod
    fail("verification", f"could not find userDefs_$USER.py in {pfw_root}")


def prompt(case: str, run_dir: Path, output_dir: Path, force: bool):
    exists = (run_dir.exists() and any(run_dir.iterdir())) or (output_dir.exists() and any(output_dir.iterdir()))
    if not exists or force:
        return
    if not sys.stdin.isatty():
        fail(case, "generated output already exists; rerun with --force")
    ans = input(f"Overwrite previous generated output for {case}? [y/N] ").strip().lower()
    if ans not in {"y", "yes"}:
        raise SystemExit("rerun cancelled")


def clean_case(run_dir: Path, output_dir: Path, force: bool):
    if run_dir.exists():
        for child in run_dir.iterdir():
            if child.is_dir() and not child.is_symlink():
                shutil.rmtree(child)
            else:
                child.unlink(missing_ok=True)
    else:
        run_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    if force:
        for child in output_dir.iterdir():
            if child.is_dir() and not child.is_symlink():
                shutil.rmtree(child)
            else:
                child.unlink(missing_ok=True)


def copy_files(pfw_root: Path, source_dir: Path, input_file: str, run_dir: Path):
    patterns = ["particleFileWriter.py", "pfw_*.py", "mpm_vv_*.py"]
    for pattern in patterns:
        for src in sorted(pfw_root.glob(pattern)):
            if src.is_file():
                shutil.copy2(src, run_dir / src.name)
    for src in sorted(source_dir.glob("*.py")):
        if src.is_file():
            shutil.copy2(src, run_dir / src.name)
    shutil.copy2(source_dir / input_file, run_dir / input_file)
    for pattern in ["*.csv", "*.txt", "*.dat", "*.png", "*.stl", "*.STL", "*.tex"]:
        for src in sorted(source_dir.glob(pattern)):
            if src.is_file():
                shutil.copy2(src, run_dir / src.name)
    copied_userdefs = []
    for src in sorted(pfw_root.glob("userDefs_*.py")):
        if src.is_file():
            dst = run_dir / src.name
            shutil.copy2(src, dst)
            copied_userdefs.append(dst)
    if copied_userdefs:
        preferred = copied_userdefs[0]
        for name in {os.environ.get("USER", ""), os.environ.get("LOGNAME", ""), getpass.getuser()}:
            if name:
                alias = run_dir / f"userDefs_{name}.py"
                if not alias.exists():
                    shutil.copy2(preferred, alias)


def append_overrides(path: Path, case: str, bank: str, walltime: str | None, submit: bool, ats_fast: bool = False):
    """Append dispatch-only controls to the staged copy of a pfw input.

    Verification source inputs own their physics, resolution, output cadence, and
    schema-visible options. The suite runner should not silently repair obsolete
    input syntax or rewrite the problem size. Use --ats-fast only for deliberately
    tiny CI reductions.
    """
    wall = f'pfw["mWallTime"] = {walltime!r}\n' if walltime else ""
    block = f"""

# -----------------------------------------------------------------------------
# Runtime dispatch controls appended by mpm_vv_case_runner.py to the staged copy
# only. Keep source verification inputs current with the GEOS/PFW schema instead
# of relying on runner compatibility rewrites.
# -----------------------------------------------------------------------------
pfw["mBatch"] = True
pfw["mBank"] = {bank!r}
{wall}pfw["mSubmitJobs"] = {bool(submit)!r}
pfw["autoRestart"] = False
"""
    if ats_fast:
        block += """
# -----------------------------------------------------------------------------
# Optional ATS-fast reduction. This path is intentionally explicit; full
# verification runs should use the resolution defined in the source input file.
# -----------------------------------------------------------------------------
def _vv_fast_int(_value, _default):
    try:
        return int(float(str(_value).strip().strip('\"').strip("'")))
    except Exception:
        return int(_default)

_vv_fast_plane = str(pfw.get("planeStrain", False)).strip().lower() in ("1", "true", "yes", "on")
_vv_fast_cpp_cap = 12 if _vv_fast_plane else 6
_vv_fast_max_partitions = 2
for _vv_key in ("xpar", "ypar", "zpar"):
    pfw[_vv_key] = max(1, min(_vv_fast_int(pfw.get(_vv_key, 1), 1), _vv_fast_max_partitions))
if _vv_fast_plane:
    pfw["zpar"] = 1

def _vv_fast_cap_cells(_nkey, _pkey, _default_cells=1):
    _p = max(1, _vv_fast_int(pfw.get(_pkey, 1), 1))
    _n = _vv_fast_int(pfw.get(_nkey, 0), 0)
    if _n <= 0:
        return max(1, _p * min(_default_cells, _vv_fast_cpp_cap))
    _cpp = max(1, (_n + _p - 1) // _p)
    return max(1, _p * min(_cpp, _vv_fast_cpp_cap))

pfw["nI"] = _vv_fast_cap_cells("nI", "xpar", _vv_fast_cpp_cap)
pfw["nJ"] = _vv_fast_cap_cells("nJ", "ypar", _vv_fast_cpp_cap)
if _vv_fast_plane:
    if "nK" in pfw:
        pfw["nK"] = max(1, min(_vv_fast_int(pfw.get("nK", 1), 1), 6))
else:
    pfw["nK"] = _vv_fast_cap_cells("nK", "zpar", 6)
pfw["mCores"] = max(1, _vv_fast_int(pfw.get("xpar", 1), 1) * _vv_fast_int(pfw.get("ypar", 1), 1) * _vv_fast_int(pfw.get("zpar", 1), 1))
if "mWallTime" not in pfw:
    pfw["mWallTime"] = "00:02:00"
"""
    text0 = path.read_text()
    if "Runtime dispatch controls appended by mpm_vv_case_runner.py" not in text0:
        path.write_text(text0.rstrip() + block + "\n")

def run_pfw(args, run_dir: Path):
    cmd = [args.python_cmd, "particleFileWriter.py", Path(args.input).name]
    log(args.case_id, "running particleFileWriter in " + str(run_dir))
    proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(proc.stdout, end="")
    (run_dir / f"{args.output_prefix}_pfw.log").write_text(proc.stdout)
    if proc.returncode != 0:
        fail(args.case_id, f"particleFileWriter failed with exit code {proc.returncode}")
    return proc.stdout


def extract_job_id(text: str) -> str:
    hits = []
    for pat in [r"Submitted job with ID\s*=\s*([0-9]+)", r"Submitted batch job\s+([0-9]+)"]:
        hits.extend(re.findall(pat, text))
    return hits[-1] if hits else ""


def submit_post(args, run_dir: Path, output_dir: Path, source_dir: Path, bank: str, geos_job: str) -> str:
    if not geos_job or shutil.which("sbatch") is None:
        return ""
    cmd = [
        args.python_cmd, "mpm_vv_postprocess.py",
        "--suite", args.suite,
        "--run-dir", str(run_dir),
        "--source-dir", str(source_dir),
        "--output-dir", str(output_dir),
        "--job-id", geos_job,
        "--case-id", args.case_id,
        "--output-prefix", args.output_prefix,
        "--python", args.python_cmd,
        "--force",
    ]
    if args.no_visit:
        cmd.append("--no-visit")
    script = run_dir / f"post_{args.output_prefix}.sbatch"
    script.write_text("\n".join([
        "#!/bin/bash",
        f"#SBATCH -A {bank}",
        f"#SBATCH -p {args.post_partition}",
        "#SBATCH -N 1",
        "#SBATCH -n 1",
        f"#SBATCH -t {args.post_walltime}",
        f"#SBATCH -J post_{args.output_prefix[:40]}",
        f"#SBATCH --dependency=afterany:{geos_job}",
        "set -euo pipefail",
        f"cd {run_dir}",
        " ".join(repr(x) for x in cmd),
        "",
    ]))
    proc = subprocess.run(["sbatch", str(script)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(proc.stdout, end="")
    m = re.search(r"Submitted batch job\s+([0-9]+)", proc.stdout)
    return m.group(1) if m else ""


def main():
    args = parse_args()
    source_dir = expand_path(args.source_dir)
    pfw_root = find_pfw_root(source_dir)
    userdefs_path, userdefs = load_userdefs(pfw_root)
    test_root = expand_path(getattr(userdefs, "testRunDirectory"))
    default_root = expand_path(getattr(userdefs, "defaultRunDirectory", test_root))
    bank = str(getattr(userdefs, "defaultBank", os.environ.get("GEOS_BANK", "")))
    if not bank:
        fail(args.case_id, "no defaultBank in userDefs and GEOS_BANK is not set")
    run_dir = test_root / args.suite / args.case_id
    output_dir = source_dir / "output" / args.case_id
    test_root.mkdir(parents=True, exist_ok=True)
    default_root.mkdir(parents=True, exist_ok=True)
    prompt(args.case_id, run_dir, output_dir, args.force)
    clean_case(run_dir, output_dir, args.force)
    copy_files(pfw_root, source_dir, args.input, run_dir)
    append_overrides(run_dir / args.input, args.case_id, bank, args.walltime, not args.no_submit, args.ats_fast)
    meta = {
        "case_id": args.case_id,
        "suite": args.suite,
        "source_dir": str(source_dir),
        "input": args.input,
        "run_dir": str(run_dir),
        "output_dir": str(output_dir),
        "userDefs": str(userdefs_path),
        "bank": bank,
    }
    run_start = time.time()
    if args.prepare_only:
        meta["status"] = "prepared"
        meta["dispatch_elapsed_seconds"] = round(time.time() - run_start, 3)
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / f"{args.output_prefix}_jobs.json").write_text(json.dumps(meta, indent=2))
        log(args.case_id, "prepared run directory: " + str(run_dir))
        return 0
    pfw_out = run_pfw(args, run_dir)
    geos_job = extract_job_id(pfw_out)
    meta["geos_job_id"] = geos_job
    if geos_job:
        log(args.case_id, "GEOS job id: " + geos_job)
    post_job = "" if args.no_post else submit_post(args, run_dir, output_dir, source_dir, bank, geos_job)
    meta["post_job_id"] = post_job
    if post_job:
        log(args.case_id, "post-process job id: " + post_job)
    meta["status"] = "submitted" if geos_job else "generated"
    meta["dispatch_elapsed_seconds"] = round(time.time() - run_start, 3)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{args.output_prefix}_jobs.json").write_text(json.dumps(meta, indent=2))
    (run_dir / f"{args.output_prefix}_jobs.json").write_text(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
