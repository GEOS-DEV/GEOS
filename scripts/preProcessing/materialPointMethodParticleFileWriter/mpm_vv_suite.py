#!/usr/bin/env python3
"""Run/report the GEOS-MPM verification suite."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


@dataclass
class Case:
    case_id: str
    source_dir: Path
    input_file: str
    output_prefix: str
    family: str
    title: str
    run_script: str = ""


def pfw_root() -> Path:
    return Path(__file__).resolve().parent


def suite_root(suite: str) -> Path:
    return pfw_root() / suite


def case_id_for(suite_dir: Path, input_path: Path) -> str:
    rel_dir = input_path.parent.relative_to(suite_dir)
    family = "_".join(rel_dir.parts) if rel_dir.parts else input_path.parent.name
    stem = input_path.stem.replace("pfw_input_", "")
    return f"{family}__{stem}"


def case_id_for_folder(suite_dir: Path, folder: Path) -> str:
    rel_dir = folder.relative_to(suite_dir)
    return "_".join(rel_dir.parts) if rel_dir.parts else folder.name


def family_label(case_id: str) -> str:
    top = case_id.split("__", 1)[0].split("_", 1)[0]
    mapping = {
        "Ftable": "F-table interpolation",
        "anneal": "Annealing",
        "bcCheck": "Boundary conditions",
        "ceramic": "Ceramic damage",
        "cohesiveZones": "Cohesive zones",
        "contact": "Contact",
        "contactBC": "Contact boundary conditions",
        "contactSurfaceGapClosure": "Contact surface/gap closure",
        "cubic": "Cubic elasticity (expected fail on this branch)",
        "diskThruPartitions": "Partition transfer",
        "expandingBar": "Expanding bar analytic comparison",
        "expandingRing": "Expanding ring analytic comparison",
        "explicitContact": "Explicit contact",
        "gas": "Gas EOS",
        "generalizedVortex": "Manufactured solution",
        "geomechanics": "Geomechanics constitutive verification",
        "implicitFluid": "Implicit fluid/pressure",
        "materialSwap": "Material swap event",
        "mmsVortex": "MMS vortex manufactured solution",
        "momentumTest": "Momentum conservation",
        "periodicBoundaries": "Periodic boundaries",
        "polymerCZ": "Polymer cohesive zone",
        "polymerHeal": "Polymer healing",
        "shrinkage": "Shrinkage",
        "sizeEffect": "Size effect",
        "spinningDisk": "Spinning disk",
        "stressControl": "Stress control",
        "temperatureTable": "Temperature table",
        "triplyPeriodic": "Triply periodic",
        "vonMisesJ": "Von Mises plasticity",
    }
    return mapping.get(top, top)


def discovery_skipped(path: Path, root: Path) -> bool:
    """Return True when a directory tree is intentionally excluded.

    Some historical verification inputs are validation-scale examples or broad
    legacy driver collections.  They remain in-tree for reference, but they
    should not be picked up by the fast verification-suite discovery pass.  A
    .mpm_vv_skip_discovery marker on any parent directory makes that policy
    explicit and local to the quarantined tree.
    """
    try:
        rel = path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    cur = root.resolve()
    if (cur / ".mpm_vv_skip_discovery").is_file():
        return True
    for part in rel.parts:
        cur = cur / part
        if cur.is_dir() and (cur / ".mpm_vv_skip_discovery").is_file():
            return True
    return False


def discover_cases(suite: str) -> list[Case]:
    root = suite_root(suite)
    out: list[Case] = []
    folder_cases: set[Path] = set()

    # New-format verification tests are folder-scoped: a single runTest script
    # owns all load cases/variants inside that folder and one postProcess*.py
    # script aggregates the folder-level results.
    for run_script in sorted(root.rglob("runTest")):
        parts = set(run_script.parts)
        if "output" in parts or "__pycache__" in parts or not run_script.is_file():
            continue
        if discovery_skipped(run_script.parent, root):
            continue
        folder = run_script.parent
        folder_cases.add(folder.resolve())
        inputs = sorted(p for p in folder.glob("pfw_input_*.py") if p.is_file())
        input_file = inputs[0].name if inputs else ""
        cid = case_id_for_folder(root, folder)
        title = cid.replace("__", " / ").replace("_", " ")
        out.append(Case(cid, folder, input_file, cid, family_label(cid), title, run_script=run_script.name))

    # Legacy fallback: until each historical directory has been migrated to the
    # folder-level contract, continue discovering standalone pfw inputs that are
    # not in a folder already managed by runTest.
    for input_path in sorted(root.rglob("pfw_input_*.py")):
        parts = set(input_path.parts)
        name = input_path.name
        if "output" in parts or "__pycache__" in parts:
            continue
        if discovery_skipped(input_path.parent, root):
            continue
        if input_path.parent.resolve() in folder_cases:
            continue
        if any(tok in name for tok in ["XXXX", "YYYY", "template", "TEMPLATE"]):
            continue
        cid = case_id_for(root, input_path)
        prefix = cid
        out.append(Case(cid, input_path.parent, input_path.name, prefix, family_label(cid), cid.replace("__", " / ").replace("_", " ")))
    return out


def select_cases(cases: list[Case], filters: list[str], max_cases: int | None) -> list[Case]:
    selected = []
    for c in cases:
        if filters:
            hay = " ".join([c.case_id, c.input_file, c.title, str(c.source_dir)])
            if not any(f in hay for f in filters):
                continue
        selected.append(c)
    if max_cases:
        selected = selected[:max_cases]
    return selected


def parse_args():
    p = argparse.ArgumentParser(description="Run/report GEOS-MPM verification suite")
    p.add_argument("action", nargs="?", default="all", choices=("all", "run", "report", "list"))
    p.add_argument("--suite", default="verification")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--max-cases", type=int)
    p.add_argument("--python", dest="python_cmd", default=os.environ.get("PFW_PYTHON", "/usr/tce/bin/python3"))
    p.add_argument("--force", "-y", action="store_true")
    p.add_argument("--stop-on-error", action="store_true")
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--no-wait", action="store_true")
    p.add_argument("--post-walltime", default='00:05:00')
    p.add_argument("--post-partition", default="pdebug")
    p.add_argument("--walltime", default=None)
    p.add_argument("--dispatch-timeout", type=float, default=None)
    p.add_argument("--jobs", type=int, default=int(os.environ.get("PFW_SUITE_JOBS", "1")), help="Number of test-folder or legacy-case dispatches to run concurrently")
    p.add_argument("--case-jobs", type=int, default=int(os.environ.get("PFW_CASE_JOBS", "1")), help="Number of subcases a folder-level runTest may dispatch concurrently")
    p.add_argument("--ats-fast", action="store_true", help="Pass explicit ultra-fast ATS sizing overrides to each staged case")
    return p.parse_args()


def run_cmd(cmd: list[str], cwd: Path, timeout=None):
    print("+ " + " ".join(cmd), flush=True)
    proc = subprocess.Popen(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    lines = []
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            lines.append(line)
        code = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        code = 124
        lines.append(f"TIMEOUT after {timeout} seconds\n")
    return code, "".join(lines)


def dispatch_one(args, runner: Path, c: Case, index: int, total: int) -> dict:
    print(f"[{args.suite}] dispatching {index}/{total}: {c.case_id}", flush=True)
    if c.run_script:
        cmd = [str(c.source_dir / c.run_script), "--suite", args.suite, "--case-id", c.case_id, "--python", args.python_cmd, "--post-walltime", args.post_walltime, "--post-partition", args.post_partition, "--jobs", str(args.case_jobs)]
    else:
        cmd = [args.python_cmd, str(runner), "--suite", args.suite, "--case-id", c.case_id, "--input", c.input_file, "--source-dir", str(c.source_dir), "--output-prefix", c.output_prefix, "--python", args.python_cmd, "--post-walltime", args.post_walltime, "--post-partition", args.post_partition]
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
    code, out = run_cmd(cmd, c.source_dir, timeout=args.dispatch_timeout)
    rec = c.__dict__.copy()
    rec["returncode"] = code
    rec["runProblem_log_excerpt"] = out[-10000:]
    outdir = c.source_dir / "output" / c.case_id
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{c.output_prefix}_runTest.log").write_text(out)
    if code != 0:
        print(f"[{args.suite}] {c.case_id} failed with {code}", flush=True)
    return rec


def dispatch(args, cases: list[Case]) -> list[dict]:
    if cases and not args.force and not args.prepare_only:
        if not sys.stdin.isatty():
            raise SystemExit("Noninteractive suite runs require --force")
        ans = input(f"Run {len(cases)} {args.suite} case(s)? [y/N] ").strip().lower()
        if ans not in {"y", "yes"}:
            raise SystemExit("cancelled")
    records = []
    runner = pfw_root() / "mpm_vv_case_runner.py"
    jobs = max(1, int(args.jobs or 1))
    if jobs == 1 or len(cases) <= 1:
        for i, c in enumerate(cases, 1):
            rec = dispatch_one(args, runner, c, i, len(cases))
            records.append(rec)
            if rec.get("returncode", 0) != 0 and args.stop_on_error:
                break
    else:
        print(f"[{args.suite}] dispatching with {jobs} concurrent worker(s)", flush=True)
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            future_to_case = {pool.submit(dispatch_one, args, runner, c, i, len(cases)): c for i, c in enumerate(cases, 1)}
            for future in as_completed(future_to_case):
                rec = future.result()
                records.append(rec)
                if rec.get("returncode", 0) != 0 and args.stop_on_error:
                    print(f"[{args.suite}] stop-on-error requested after {rec.get('case_id')}", flush=True)
    report_dir = suite_root(args.suite) / f"{args.suite}_suite_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / f"{args.suite}_suite_dispatch_records.json").write_text(json.dumps(records, indent=2, default=str))
    return records

def wait_for_posts(args, cases: list[Case]) -> None:
    if args.no_wait or shutil.which("squeue") is None:
        return
    job_ids = []
    for c in cases:
        for path in (c.source_dir / "output" / c.case_id).glob("*_jobs.json"):
            try:
                data = json.loads(path.read_text())
            except Exception:
                continue
            jid = str(data.get("post_job_id") or "")
            if jid:
                job_ids.append(jid)
    if not job_ids:
        return
    print(f"[{args.suite}] waiting for {len(job_ids)} post-process job(s)", flush=True)
    import time
    remaining = set(job_ids)
    while remaining:
        proc = subprocess.run(["squeue", "-h", "-j", ",".join(sorted(remaining)), "-o", "%i"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        running = set(x.strip().split("_")[0] for x in proc.stdout.splitlines() if x.strip())
        remaining &= running
        if remaining:
            print(f"[{args.suite}] still waiting on {len(remaining)} job(s): " + ",".join(sorted(remaining)), flush=True)
            time.sleep(20)


def sacct_row(job_id: str) -> dict:
    if not job_id or shutil.which("sacct") is None:
        return {}
    proc = subprocess.run(["sacct", "-j", job_id, "-X", "-n", "-P", "-o", "JobID,State,ExitCode,Elapsed,TotalCPU,CPUTimeRAW,NCPUS"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    line = next((l for l in proc.stdout.splitlines() if l.strip()), "")
    if not line:
        return {}
    cols = line.split("|")
    keys = ["JobID", "State", "ExitCode", "Elapsed", "TotalCPU", "CPUTimeRAW", "NCPUS"]
    return {k: cols[i] if i < len(cols) else "" for i, k in enumerate(keys)}


def split_job_ids(value: object) -> list[str]:
    """Extract scheduler job ids from a scalar or comma-separated metadata field."""
    if value is None:
        return []
    return re.findall(r"\d+", str(value))


def scheduler_failure_state(state: str) -> bool:
    return str(state or "").startswith(("FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL", "OUT_OF_MEMORY"))


def timing_rows_for_status(data: dict) -> list[dict]:
    """Build one timing row per concrete GEOS subjob.

    Folder-level tests, such as mmsVortex, have several subcases under one
    report section. Legacy tests usually have a single GEOS job directly on the
    case record. In both forms, GEOS wall/CPU data come from sacct when
    available; dispatch_elapsed_seconds is local runner wall time for staging,
    PFW generation, and scheduler submission.
    """
    rows: list[dict] = []
    subcases = data.get("subcases", []) or []
    if subcases:
        for sub in subcases:
            jobs = sub.get("jobs", {}) or {}
            job_ids = split_job_ids(jobs.get("geos_job_id", ""))
            geos_job_id = job_ids[0] if job_ids else str(jobs.get("geos_job_id", "") or "")
            geos = sacct_row(geos_job_id) if geos_job_id else {}
            rows.append(
                {
                    "case_id": data.get("case_id", ""),
                    "subcase": sub.get("label") or sub.get("name") or sub.get("case_id", ""),
                    "geos_job_id": geos_job_id,
                    "geos_state": geos.get("State", ""),
                    "geos_elapsed": geos.get("Elapsed", ""),
                    "geos_total_cpu": geos.get("TotalCPU", ""),
                    "geos_ncpus": geos.get("NCPUS", ""),
                    "dispatch_elapsed_seconds": jobs.get("dispatch_elapsed_seconds", sub.get("dispatch_elapsed_seconds", "")),
                    "returncode": sub.get("returncode", ""),
                }
            )
        return rows

    job_ids = split_job_ids(data.get("geos_job_id", ""))
    geos_job_id = job_ids[0] if job_ids else str(data.get("geos_job_id", "") or "")
    geos = sacct_row(geos_job_id) if geos_job_id else {}
    rows.append(
        {
            "case_id": data.get("case_id", ""),
            "subcase": "",
            "geos_job_id": geos_job_id,
            "geos_state": geos.get("State", ""),
            "geos_elapsed": geos.get("Elapsed", ""),
            "geos_total_cpu": geos.get("TotalCPU", ""),
            "geos_ncpus": geos.get("NCPUS", ""),
            "dispatch_elapsed_seconds": data.get("dispatch_elapsed_seconds", ""),
            "returncode": data.get("returncode", ""),
        }
    )
    return rows


def case_job_summary(status: dict) -> str:
    rows = status.get("timing_rows", []) or []
    ids = [str(row.get("geos_job_id", "")) for row in rows if row.get("geos_job_id")]
    if ids:
        return ", ".join(ids)
    return str(status.get("geos_job_id", "") or "")


def case_state_summary(status: dict) -> str:
    rows = status.get("timing_rows", []) or []
    states = [str(row.get("geos_state", "")) for row in rows if row.get("geos_state")]
    if states:
        unique = []
        for state in states:
            if state not in unique:
                unique.append(state)
        return ", ".join(unique)
    return str((status.get("geos_sacct", {}) or {}).get("State", ""))


def case_elapsed_summary(status: dict) -> str:
    rows = status.get("timing_rows", []) or []
    elapsed = []
    for row in rows:
        wall = str(row.get("geos_elapsed", "") or "")
        if wall:
            label = str(row.get("subcase", "") or row.get("case_id", ""))
            elapsed.append(f"{label}: {wall}" if label else wall)
    if elapsed:
        return "; ".join(elapsed)
    dispatch = [str(row.get("dispatch_elapsed_seconds", "")) for row in rows if row.get("dispatch_elapsed_seconds") not in (None, "")]
    if dispatch:
        return "dispatch " + ", ".join(dispatch) + " s"
    return ""






def collect_case_status(c: Case) -> dict:
    outdir = c.source_dir / "output" / c.case_id
    data = c.__dict__.copy()
    data.update({"output_dir": str(outdir), "status": "not run", "issue": "", "warnings": []})

    # Recover the most recent job metadata if report generation is run without a
    # fresh dispatch.  This also prevents stale dispatch-failure logs from
    # overriding later real runs.
    for jobs in sorted(outdir.glob("*_jobs.json")):
        try:
            data.update(json.loads(jobs.read_text()))
        except Exception:
            pass

    geos_job_ids = split_job_ids(data.get("geos_job_id", ""))
    geos = sacct_row(geos_job_ids[0]) if len(geos_job_ids) == 1 else {}
    post_job_ids = split_job_ids(data.get("post_job_id", ""))
    post = sacct_row(post_job_ids[0]) if len(post_job_ids) == 1 else {}
    data["geos_sacct"] = geos
    data["post_sacct"] = post
    data["timing_rows"] = timing_rows_for_status(data)

    logs = []
    for p in sorted(outdir.glob("*.log")):
        try:
            logs.append((p.name, p.read_text(errors="replace")[-25000:]))
        except Exception:
            pass
    data["logs"] = [{"name": n, "text": t} for n, t in logs]

    pngs = list(outdir.glob("*.png"))
    csvs = list(outdir.glob("*.csv"))
    json_products = [p for p in outdir.glob("*.json") if not p.name.endswith("_jobs.json")]
    tex_products = list(outdir.glob("*_results.tex"))
    timing_rows = data.get("timing_rows", []) or []
    completed_timing_rows = [row for row in timing_rows if row.get("geos_state") == "COMPLETED"]
    geos_job_evidence = bool(data.get("geos_job_id") or completed_timing_rows or any(row.get("geos_job_id") for row in timing_rows))
    job_or_product_evidence = bool(pngs or csvs or json_products or tex_products or geos_job_evidence or data.get("post_job_id"))
    subcases = data.get("subcases", []) or []
    subcase_failures = [s for s in subcases if int(s.get("returncode", 0) or 0) != 0]
    subcase_scheduler_failures = [row for row in data.get("timing_rows", []) if scheduler_failure_state(row.get("geos_state", ""))]

    hard_log_names = (
        "geos_slurm_check", "pfw", "runProblem", "particleFileWriter",
        "slurm", "visit_render", "postprocess", "postProcess",
    )
    soft_plot_log = re.compile(r"plot(Reactions|Box|History|.*)?|legacy", re.I)

    hard_text_parts = []
    soft_text_parts = []
    for name, text in logs:
        if "Unknown case" in text and job_or_product_evidence:
            # Stale logs from the early suite-dispatch bug; keep them in the
            # report diagnostics but do not let them determine case status.
            data["warnings"].append(f"ignored stale dispatch log: {name}")
            continue
        if any(k in name for k in hard_log_names) and not soft_plot_log.search(name):
            hard_text_parts.append(text)
        elif "geos_slurm_check" in name or "runProblem" in name or "pfw" in name:
            hard_text_parts.append(text)
        else:
            soft_text_parts.append(text)
    hard_text = "\n".join(hard_text_parts)
    soft_text = "\n".join(soft_text_parts)
    all_text = "\n".join(t for _, t in logs)

    expected_re = re.compile(r"ElasticCubic.*invalid within Constitutive|tag \"?ElasticCubic\"? is invalid", re.I)
    hard_fail_re = re.compile(
        r"XML parsing error|EXCEPTION|contains unused attribute|missing required attribute|invalid within|"
        r"Floating point error|TIME LIMIT|NODE_FAIL|OUT_OF_MEMORY|particleFileWriter failed|"
        r"\bfailure_marker\b|\bERROR\b",
        re.I,
    )
    cancelled_re = re.compile(r"CANCELLED|TIMEOUT|DUE to SIGNAL Killed", re.I)
    soft_fail_re = re.compile(r"Traceback|legacy plot failed|plot script failed", re.I)

    if expected_re.search(all_text):
        data["status"] = "expected-fail"
        data["issue"] = "ElasticCubic is not registered in this branch"
    elif subcase_failures:
        data["status"] = "failed"
        data["issue"] = f"{len(subcase_failures)} subcase dispatch failure(s)"
    elif subcase_scheduler_failures:
        data["status"] = "failed"
        data["issue"] = ", ".join(
            f"{row.get('subcase') or row.get('case_id')}: {row.get('geos_state')}"
            for row in subcase_scheduler_failures[:3]
        )
    elif hard_fail_re.search(hard_text) or cancelled_re.search(hard_text):
        data["status"] = "failed"
        m = hard_fail_re.search(hard_text) or cancelled_re.search(hard_text)
        data["issue"] = m.group(0) if m else "failure marker"
    elif geos.get("State", "").startswith(("FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL", "OUT_OF_MEMORY")):
        data["status"] = "failed"
        data["issue"] = geos.get("State", "scheduler failure")
    elif pngs:
        if soft_fail_re.search(soft_text):
            data["status"] = "passed-with-warnings"
            data["issue"] = "optional legacy plotting script failed"
            data["warnings"].append("Generic products exist; optional legacy plotting traceback ignored for pass/fail status.")
        else:
            data["status"] = "passed"
    elif pngs or csvs or json_products or tex_products:
        if soft_fail_re.search(soft_text):
            data["status"] = "passed-with-warnings"
            data["issue"] = "verification products exist; optional plotting warning found"
            data["warnings"].append("Verification products exist; optional plotting traceback ignored for pass/fail status.")
        else:
            data["status"] = "passed"
    elif geos.get("State") == "COMPLETED" or geos_job_evidence:
        if soft_fail_re.search(soft_text):
            data["status"] = "completed-no-figures"
            data["issue"] = "GEOS completed; optional plot script failed and no PNG products were found"
        else:
            data["status"] = "completed-no-figures"
            data["issue"] = "GEOS completed but no PNG products were found"
    elif soft_fail_re.search(soft_text):
        data["status"] = "postprocess-warning"
        data["issue"] = "optional legacy plotting script failed"
    else:
        data["status"] = "not run"
    return data

def latex_escape(text: object) -> str:
    s = str(text if text is not None else "")
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(repl.get(ch, ch) for ch in s)


def include_graphic(path: Path, width="0.31\\linewidth") -> str:
    return r"\includegraphics[width=" + width + "]{" + latex_escape(path.as_posix()) + "}"


def write_schematic(c: Case, report_dir: Path) -> Path:
    out = report_dir / "schematics" / f"{c.case_id}_schematic.png"
    script = pfw_root() / "mpm_vv_schematic.py"
    if not out.is_file():
        subprocess.run([sys.executable, str(script), "--input", str(c.source_dir / c.input_file), "--output", str(out), "--case-id", c.case_id, "--pfw-root", str(pfw_root())], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out


def case_description(c: Case) -> str:
    fam = c.family
    title = c.title
    return f"This verification problem exercises {fam.lower()} using the mesh, particle seeding, material model, boundary conditions, and history outputs defined in {c.input_file}. The plots below include initial and final state renderings plus the reaction and box-average histories needed for analytical comparison."


def features_for(c: Case) -> list[str]:
    text = (c.source_dir / c.input_file).read_text(errors="replace") if (c.source_dir / c.input_file).is_file() else ""
    features = []
    for label, pat in [
        ("FMPM update", "FMPM"),
        ("FLIP/PIC/XPIC transfer", "FLIP|PIC|XPIC"),
        ("periodic boundaries", "periodic"),
        ("reaction-history output", "reactionHistory"),
        ("box-average history", "boxAverageHistory"),
        ("cohesive-zone fields", "CohesiveZone|cohesive"),
        ("Geomechanics constitutive model", "Geomechanics|ghareb"),
        ("ElasticIsotropic material", "ElasticIsotropic|elastic"),
        ("VonMisesJ plasticity", "VonMisesJ"),
        ("temperature or annealing event", "TemperatureProfile|Anneal|temperature"),
        ("contact/multi-field interaction", "Contact|contact|ContactGroup"),
    ]:
        if re.search(pat, text, re.I):
            features.append(label)
    return features or ["MPM particle generation and explicit dynamics"]


def plot_plan_latex(outdir: Path, prefix: str) -> str:
    plan_path = outdir / f"{prefix}_plot_requirements.json"
    if not plan_path.is_file():
        return "No plot plan JSON was generated."
    try:
        plan = json.loads(plan_path.read_text())
    except Exception:
        return "Plot plan JSON could not be read."
    rows = []
    for p in plan.get("plots", []):
        rows.append(" & ".join(latex_escape("; ".join(p.get(k, [])) if isinstance(p.get(k), list) else p.get(k, "")) for k in ["script", "data_sources", "variables", "processing"]) + r" \\")
    return r"""
{\scriptsize
\begin{tabularx}{\linewidth}{p{0.18\linewidth}X X X}
\toprule
Script/source & Data source & Variables & Processing \\
\midrule
""" + "\n".join(rows[:12]) + "\n\\bottomrule\n\\end{tabularx}\n}\n"


def diagnostics_latex(status: dict, max_chars=5000) -> str:
    logs = status.get("logs", [])
    important = []
    pat = re.compile(r"Unknown case|XML parsing error|EXCEPTION|Error cause|missing required attribute|contains unused attribute|invalid within|Floating point error|Traceback|CANCELLED|Killed|TIME LIMIT", re.I)
    for log in logs:
        lines = [line for line in log.get("text", "").splitlines() if pat.search(line)]
        if lines:
            important.append(f"--- {log.get('name')} ---\n" + "\n".join(lines[-60:]))
    if not important:
        return ""
    text = "\n".join(important)[-max_chars:]
    return r"\paragraph{Diagnostics.}\begin{lstlisting}[basicstyle=\tiny\ttfamily,breaklines=true]" + "\n" + text.replace("\x0c", " ") + "\n" + r"\end{lstlisting}"


def write_report(args, cases: list[Case]) -> None:
    report_dir = suite_root(args.suite) / f"{args.suite}_suite_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    statuses = [collect_case_status(c) for c in cases]
    (report_dir / f"{args.suite}_suite_status.json").write_text(json.dumps(statuses, indent=2, default=str))
    tex = []
    tex.append(r"\documentclass[10pt]{article}")
    tex.append(r"\usepackage[margin=0.55in]{geometry}")
    tex.append(r"\usepackage{graphicx,booktabs,longtable,tabularx,array,pdflscape,xcolor,hyperref,listings,amsmath}")
    tex.append(r"\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}")
    tex.append(r"\setlength{\parindent}{0pt}")
    tex.append(r"\begin{document}")
    tex.append(r"\title{GEOS-MPM Verification Suite Report}")
    tex.append(r"\author{Generated by verification/runSuite}")
    tex.append(r"\date{" + latex_escape(datetime.now().strftime("%B %d, %Y")) + "}")
    tex.append(r"\maketitle\tableofcontents\clearpage")
    tex.append(r"\section{Running and rerunning}")
    tex.append(r"Cases are staged under \texttt{testRunDirectory/verification/<case>}. The suite records failures and continues unless \texttt{--stop-on-error} is supplied. Use \texttt{--force} for noninteractive reruns.")
    tex.append(r"\begin{lstlisting}[basicstyle=\small\ttfamily]" + "\ncd scripts/preProcessing/materialPointMethodParticleFileWriter\nverification/runSuite run --python /usr/tce/bin/python3 --force --jobs 4 --case-jobs 2\nverification/runSuite report --python /usr/tce/bin/python3\n" + r"\end{lstlisting}")
    tex.append(r"\section{Timing and dispatch summary}")
    tex.append(r"\begin{landscape}{\tiny\begin{longtable}{p{0.25\linewidth}p{0.10\linewidth}p{0.19\linewidth}p{0.15\linewidth}p{0.10\linewidth}p{0.11\linewidth}p{0.13\linewidth}}\toprule Case & Status & Family & GEOS job(s) & GEOS state & GEOS wall & Issue \\ \midrule\endhead")
    for st in statuses:
        tex.append(" & ".join(latex_escape(x) for x in [st.get("case_id", ""), st.get("status", ""), st.get("family", ""), case_job_summary(st), case_state_summary(st), case_elapsed_summary(st), st.get("issue", "")]) + r" \\")
    tex.append(r"\bottomrule\end{longtable}}\end{landscape}")

    tex.append(r"\paragraph{Subcase timing detail.} Dispatch time is local wall time spent staging inputs, running PFW, and submitting the GEOS job. GEOS wall and CPU columns are populated from \texttt{sacct} when available.")
    tex.append(r"\begin{landscape}{\tiny\begin{longtable}{p{0.24\linewidth}p{0.18\linewidth}p{0.10\linewidth}p{0.10\linewidth}p{0.10\linewidth}p{0.10\linewidth}p{0.07\linewidth}p{0.09\linewidth}}\toprule Case & Subcase & GEOS job & State & GEOS wall & Total CPU & NCPUS & Dispatch (s) \\ \midrule\endhead")
    for st in statuses:
        for row in st.get("timing_rows", []) or []:
            tex.append(" & ".join(latex_escape(x) for x in [row.get("case_id", ""), row.get("subcase", ""), row.get("geos_job_id", ""), row.get("geos_state", ""), row.get("geos_elapsed", ""), row.get("geos_total_cpu", ""), row.get("geos_ncpus", ""), row.get("dispatch_elapsed_seconds", "")]) + r" \\")
    tex.append(r"\bottomrule\end{longtable}}\end{landscape}")
    failures = [s for s in statuses if s.get("status") in {"failed", "expected-fail", "completed-no-figures"}]
    if failures:
        tex.append(r"\section{Failures, expected failures, and diagnostics}")
        for st in failures:
            tex.append(r"\subsection{" + latex_escape(st.get("case_id", "")) + "}")
            tex.append(r"\textbf{Status:} " + latex_escape(st.get("status", "")) + r"\quad \textbf{Issue:} " + latex_escape(st.get("issue", "")))
            tex.append(diagnostics_latex(st, max_chars=3500))
    tex.append(r"\section{Case studies}")
    for c, st in zip(cases, statuses):
        outdir = c.source_dir / "output" / c.case_id
        prefix = c.output_prefix
        tex.append(r"\subsection{" + latex_escape(c.title) + "}")
        case_tex = c.source_dir / "test.tex"
        if case_tex.is_file():
            rel_case_tex = Path(os.path.relpath(case_tex, report_dir))
            rel_source_dir = Path(os.path.relpath(c.source_dir, report_dir))
            rel_output_dir = Path(os.path.relpath(outdir, report_dir))
            tex.append(r"\begingroup")
            tex.append(r"\def\CaseId{" + latex_escape(c.case_id) + r"}")
            tex.append(r"\def\CaseSourceDir{" + str(rel_source_dir).replace("\\", "/") + r"}")
            tex.append(r"\def\CaseOutputDir{" + str(rel_output_dir).replace("\\", "/") + r"}")
            tex.append(r"\input{" + str(rel_case_tex).replace("\\", "/") + r"}")
            case_tex_text = case_tex.read_text(errors="replace")
            if "_results.tex" not in case_tex_text:
                # Split-folder placeholders often rely on the common generic
                # post-processor.  Include its generated fragment here so every
                # current-standard folder can show at least one history plot and
                # one VisIt smoke-render when products are available.
                generic_results = sorted(outdir.glob("*_results.tex"))
                if generic_results:
                    rel_result = Path(os.path.relpath(generic_results[0], report_dir))
                    tex.append(r"\input{" + str(rel_result).replace("\\", "/") + r"}")
            tex.append(r"\endgroup")
            if st.get("status") in {"failed", "expected-fail", "completed-no-figures"}:
                tex.append(diagnostics_latex(st, max_chars=3500))
            tex.append(r"\clearpage")
            continue
        tex.append(latex_escape(case_description(c)))
        schematic = write_schematic(c, report_dir)
        if schematic.is_file():
            tex.append(r"\begin{figure}[h!]\centering " + include_graphic(schematic.relative_to(report_dir), "0.70\\linewidth") + r"\caption{Geometry-aware schematic for " + latex_escape(c.case_id) + r".}\end{figure}")
        tex.append(r"\paragraph{Code features exercised.}")
        tex.append(r"\begin{itemize}")
        for f in features_for(c):
            tex.append(r"\item " + latex_escape(f))
        tex.append(r"\end{itemize}")
        tex.append(r"\paragraph{Verification plotting plan.}")
        plan_path = outdir / f"{prefix}_plot_requirements.json"
        if not plan_path.is_file():
            try:
                sys.path.insert(0, str(pfw_root()))
                import mpm_vv_plot_tools as _vv_tools
                _vv_tools.write_plan(plan_path, _vv_tools.make_plot_plan(c.source_dir, c.case_id, c.family))
            except Exception:
                pass
        tex.append(plot_plan_latex(outdir, prefix))
        pngs = sorted(outdir.glob("*.png"))
        # Prefer initial/final/history and then any additional verification curves.
        def score(p: Path):
            n = p.name.lower()
            if "initial" in n: return (0, n)
            if "final" in n: return (1, n)
            if "reaction" in n or "history" in n: return (2, n)
            if "box" in n or "stress" in n or "pressure" in n or "strain" in n: return (3, n)
            return (4, n)
        pngs = sorted(pngs, key=score)
        if pngs:
            tex.append(r"\begin{figure}[h!]\centering")
            for p in pngs[:6]:
                tex.append(include_graphic(p.relative_to(report_dir), "0.31\\linewidth") if False else include_graphic(Path(os.path.relpath(p, report_dir)), "0.31\\linewidth"))
            tex.append(r"\caption{Initial/final state and verification history plots for " + latex_escape(c.case_id) + r".}\end{figure}")
        else:
            tex.append(r"\begin{center}\fbox{missing initial/final/history plots}\end{center}")
        if st.get("status") in {"failed", "expected-fail", "completed-no-figures"}:
            tex.append(diagnostics_latex(st, max_chars=3500))
        tex.append(r"\clearpage")
    tex.append(r"\end{document}")
    tex_path = report_dir / f"{args.suite}_suite_report.tex"
    tex_path.write_text("\n".join(tex))
    for _ in range(2):
        subprocess.run(["pdflatex", "-interaction=nonstopmode", tex_path.name], cwd=report_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(f"[{args.suite}] report written to {report_dir / (args.suite + '_suite_report.pdf')}")


def main():
    args = parse_args()
    cases = select_cases(discover_cases(args.suite), args.case, args.max_cases)
    print(f"[{args.suite}] {len(cases)} selected case(s)")
    if args.action == "list":
        for c in cases:
            print(f"{c.case_id}\t{c.input_file}\t{c.source_dir}\t{c.run_script or 'legacy'}")
        return 0
    if args.action in {"all", "run"}:
        dispatch(args, cases)
        if not args.prepare_only:
            wait_for_posts(args, cases)
    if args.action in {"all", "report"}:
        write_report(args, cases)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
