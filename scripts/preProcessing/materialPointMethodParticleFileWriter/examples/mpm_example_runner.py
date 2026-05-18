#!/usr/bin/env python3
"""Stage and submit one MPM example case from an examples/<case>/runProblem wrapper."""
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
from datetime import datetime
from pathlib import Path


def log(case: str, msg: str) -> None:
    print(f"[{case}] {msg}")


def fail(case: str, msg: str) -> None:
    print(f"[{case}] ERROR: {msg}", file=sys.stderr)
    raise SystemExit(1)


def expand_path(path: str | os.PathLike) -> Path:
    return Path(os.path.expandvars(os.path.expanduser(str(path)))).resolve()


def parse_args():
    p = argparse.ArgumentParser(description="Stage, generate, and submit one PFW example")
    p.add_argument("--case-id", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--source-dir", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--initial-variable", default="Damage")
    p.add_argument("--final-variable", default="Damage")
    p.add_argument("--initial-range-mode", choices=("unit", "auto"), default="auto")
    p.add_argument("--final-range-mode", choices=("unit", "auto"), default="auto")
    p.add_argument("--view", default="xy")
    p.add_argument("--mesh", default="CellRegion1")
    p.add_argument("--python", dest="python_cmd", default=os.environ.get("PFW_PYTHON", "/usr/tce/bin/python3"))
    p.add_argument("--force", "-y", action="store_true")
    p.add_argument("--prepare-only", action="store_true", help="stage files and write copied input, but do not run PFW")
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--post-walltime", default="00:05:00")
    p.add_argument("--post-partition", default="pdebug")
    p.add_argument("--walltime", default=None, help="Optional GEOS walltime override; by default use pfw_input mWallTime")
    p.add_argument("--no-submit", action="store_true", help="run PFW without submitting GEOS")
    return p.parse_args()


def find_pfw_root(source_dir: Path) -> Path:
    for parent in [source_dir] + list(source_dir.parents):
        if (parent / "particleFileWriter.py").is_file() and (parent / "pfw_geometryObjects.py").is_file():
            return parent
    fail("mpm_example", f"could not locate PFW root above {source_dir}")


def load_userdefs(pfw_root: Path):
    candidates = []
    user = os.environ.get("USER", "")
    if user:
        candidates.append(pfw_root / f"userDefs_{user}.py")
    candidates += sorted(pfw_root.glob("userDefs_*.py"))
    for path in candidates:
        if path.is_file():
            spec = importlib.util.spec_from_file_location("_pfw_userdefs", path)
            module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(module)
            return path, module
    fail("mpm_example", f"could not find userDefs_$USER.py in {pfw_root}")


def prompt_to_overwrite(case: str, run_dir: Path, output_dir: Path, prefix: str, force: bool) -> None:
    generated = []
    if run_dir.exists() and any(run_dir.iterdir()):
        generated.append(str(run_dir))
    if output_dir.exists():
        generated += [str(p) for p in output_dir.glob(prefix + "*")]
    if not generated:
        return
    if force:
        return
    if not sys.stdin.isatty():
        fail(case, "generated output already exists; rerun with --force in noninteractive mode")
    print(f"Existing generated output for {case}:")
    for item in generated[:20]:
        print(f"  {item}")
    if len(generated) > 20:
        print(f"  ... and {len(generated)-20} more")
    ans = input(f"Overwrite previous {case} generated output and rerun? [y/N] ").strip().lower()
    if ans not in {"y", "yes"}:
        raise SystemExit("rerun cancelled")


def clean_case_outputs(run_dir: Path, output_dir: Path, prefix: str, force: bool) -> None:
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
        for pattern in [prefix + "*.png", prefix + "*.csv", prefix + "*.log", prefix + "*.json", prefix + "*.txt"]:
            for path in output_dir.glob(pattern):
                if path.is_file() or path.is_symlink():
                    path.unlink(missing_ok=True)
                elif path.is_dir():
                    shutil.rmtree(path)


def copy_pfw_files(pfw_root: Path, run_dir: Path, source_dir: Path, input_file: str) -> None:
    patterns = ["particleFileWriter.py", "pfw_*.py", "mpm_example_postprocess.py"]
    copied = set()
    for pattern in patterns:
        for src in pfw_root.glob(pattern):
            if src.is_file():
                shutil.copy2(src, run_dir / src.name)
                copied.add(src.name)
    examples_postprocess = pfw_root / "examples" / "mpm_example_postprocess.py"
    if examples_postprocess.is_file():
        shutil.copy2(examples_postprocess, run_dir / examples_postprocess.name)
    for src in source_dir.glob("*.py"):
        if src.is_file() and src.name != input_file:
            shutil.copy2(src, run_dir / src.name)
    shutil.copy2(source_dir / input_file, run_dir / input_file)
    # Copy obvious non-Python dependencies beside the input.  PFW dependency tags still handle
    # cross-directory dependencies such as STL files or pfw_materials.py.
    for pattern in ["*.stl", "*.STL", "*.csv", "*.txt"]:
        for src in source_dir.glob(pattern):
            if src.is_file():
                shutil.copy2(src, run_dir / src.name)

    # The staged particleFileWriter.py imports userDefs_<user> from the run directory.
    # Copy the local userDefs files so staged runs do not depend on running from the
    # source tree and so getpass/USER differences in batch mode are tolerated.
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


def append_runtime_overrides(input_path: Path, case: str, bank: str, geos_path: str, walltime: str, submit: bool) -> None:
    walltime_line = f'pfw["mWallTime"] = {walltime!r}\n' if walltime else ""
    block = f"""

# -----------------------------------------------------------------------------
# Runtime overrides appended by examples/mpm_example_runner.py in the copied
# Lustre run directory only.  The source example input remains copy-and-modify.
# -----------------------------------------------------------------------------
pfw["mBatch"] = True
pfw["mBank"] = {bank!r}
{walltime_line}pfw["mSubmitJobs"] = {bool(submit)!r}
pfw["autoRestart"] = False
pfw["outputType"] = "silo"
try:
    _pfw_end_time = float(pfw.get("endTime", stopTime))
except Exception:
    _pfw_end_time = 1.0
pfw["plotInterval"] = _pfw_end_time
pfw["restartInterval"] = 2.0 * _pfw_end_time
_interp = pfw.get("fTableInterpType", "Linear")
if _interp == 0 or _interp == "0":
    pfw["fTableInterpType"] = "Linear"
elif _interp == 1 or _interp == "1":
    pfw["fTableInterpType"] = "Cosine"
elif _interp == 2 or _interp == "2":
    pfw["fTableInterpType"] = "Smoothstep"
elif _interp not in ("Linear", "Cosine", "Smoothstep"):
    pfw["fTableInterpType"] = "Linear"
_contact_gap = pfw.get("contactGapCorrection", None)
if _contact_gap == 0 or _contact_gap == "0":
    pfw["contactGapCorrection"] = "Simple"
elif _contact_gap == 1 or _contact_gap == "1":
    pfw["contactGapCorrection"] = "Implicit"
elif _contact_gap == 2 or _contact_gap == "2":
    pfw["contactGapCorrection"] = "Softened"
elif isinstance(_contact_gap, str) and _contact_gap.lower() in ("simple", "implicit", "softened"):
    pfw["contactGapCorrection"] = _contact_gap[:1].upper() + _contact_gap[1:].lower()
"""
    text = input_path.read_text()
    if "Runtime overrides appended by examples/mpm_example_runner.py" not in text:
        input_path.write_text(text.rstrip() + block + "\n")


def run_pfw(args, run_dir: Path):
    input_stem = Path(args.input).stem
    cmd = [args.python_cmd, "particleFileWriter.py", input_stem]
    log(args.case_id, "running particleFileWriter in " + str(run_dir))
    proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(proc.stdout, end="")
    (run_dir / f"{args.output_prefix}_pfw.log").write_text(proc.stdout)
    if proc.returncode != 0:
        fail(args.case_id, f"particleFileWriter failed with exit code {proc.returncode}")
    return proc.stdout


def extract_geos_job_id(text: str) -> str:
    patterns = [r"Submitted job with ID\s*=\s*([0-9]+)", r"Submitted batch job\s+([0-9]+)", r"GEOS job id:\s*([0-9]+)"]
    found = []
    for pat in patterns:
        found += re.findall(pat, text)
    return found[-1] if found else ""


def submit_postprocess(args, run_dir: Path, output_dir: Path, bank: str, geos_job: str) -> str:
    if not geos_job:
        log(args.case_id, "no GEOS job id detected; skipping dependent post-process submission")
        return ""
    sbatch = shutil.which("sbatch")
    if sbatch is None:
        log(args.case_id, "sbatch not found; skipping post-process submission")
        return ""
    cmd = [
        args.python_cmd,
        "mpm_example_postprocess.py",
        "--run-dir", str(run_dir),
        "--output-dir", str(output_dir),
        "--job-id", geos_job,
        "--case-name", args.case_id,
        "--output-prefix", args.output_prefix,
        "--initial-variable", args.initial_variable,
        "--final-variable", args.final_variable,
        "--initial-range-mode", args.initial_range_mode,
        "--final-range-mode", args.final_range_mode,
        "--view", args.view,
        "--mesh", args.mesh,
        "--force",
    ]
    if args.no_visit:
        cmd.append("--no-visit")
    wrap = "cd " + str(run_dir) + " && " + " ".join(shlex_quote(c) for c in cmd)
    sbatch_cmd = [
        sbatch,
        "-A", bank,
        "-p", args.post_partition,
        "-N", "1",
        "-n", "1",
        "-t", args.post_walltime,
        "-J", "post_" + args.case_id[:40],
        "--dependency=afterany:" + geos_job,
        "--wrap", wrap,
    ]
    proc = subprocess.run(sbatch_cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(proc.stdout, end="")
    if proc.returncode != 0:
        fail(args.case_id, "post-process sbatch failed")
    m = re.findall(r"Submitted batch job\s+([0-9]+)", proc.stdout)
    return m[-1] if m else ""


def shlex_quote(s: str) -> str:
    import shlex
    return shlex.quote(str(s))


def main():
    args = parse_args()
    source_dir = expand_path(args.source_dir)
    pfw_root = find_pfw_root(source_dir)
    userdefs_path, userdefs = load_userdefs(pfw_root)
    test_dir = expand_path(getattr(userdefs, "testRunDirectory", f"/p/lustre1/{os.environ.get('USER','user')}/geosxTests"))
    default_dir = expand_path(getattr(userdefs, "defaultRunDirectory", f"/p/lustre1/{os.environ.get('USER','user')}/geosxRuns"))
    bank = str(getattr(userdefs, "defaultBank", os.environ.get("GEOS_BANK", ""))).strip()
    if not bank:
        fail(args.case_id, "defaultBank is empty; update userDefs_$USER.py or pass GEOS_BANK")
    geos_path = str(expand_path(getattr(userdefs, "geosPath", "")))

    run_dir = test_dir / args.case_id
    output_dir = source_dir / "output"
    test_dir.mkdir(parents=True, exist_ok=True)
    default_dir.mkdir(parents=True, exist_ok=True)

    log(args.case_id, "PFW root: " + str(pfw_root))
    log(args.case_id, "userDefs: " + str(userdefs_path))
    log(args.case_id, "python: " + str(args.python_cmd))
    log(args.case_id, "GEOS executable from userDefs: " + geos_path)
    log(args.case_id, "run directory: " + str(run_dir))
    log(args.case_id, "example output directory: " + str(output_dir))
    log(args.case_id, "default bank/account: " + bank)

    prompt_to_overwrite(args.case_id, run_dir, output_dir, args.output_prefix, args.force)
    clean_case_outputs(run_dir, output_dir, args.output_prefix, args.force)
    copy_pfw_files(pfw_root, run_dir, source_dir, args.input)
    append_runtime_overrides(run_dir / args.input, args.case_id, bank, geos_path, args.walltime, not args.no_submit)

    if args.prepare_only:
        log(args.case_id, "prepare-only requested; staged files and copied input are ready")
        return

    pfw_stdout = run_pfw(args, run_dir)
    geos_job = extract_geos_job_id(pfw_stdout)
    if geos_job:
        log(args.case_id, "GEOS job id: " + geos_job)
    else:
        log(args.case_id, "GEOS job id: not detected")
    post_job = submit_postprocess(args, run_dir, output_dir, bank, geos_job)
    if post_job:
        log(args.case_id, "post-process job id: " + post_job)

    payload = {
        "case": args.case_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "jobs": {"geos": geos_job, "post": post_job},
        "run_dir": str(run_dir),
        "output_dir": str(output_dir),
        "input": args.input,
    }
    for path in [run_dir / f"{args.output_prefix}_jobs.json", output_dir / f"{args.output_prefix}_jobs.json"]:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
