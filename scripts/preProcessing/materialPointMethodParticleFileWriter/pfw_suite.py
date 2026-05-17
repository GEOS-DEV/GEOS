#!/usr/bin/env python3
"""
Suite runner and reporter for GEOS-MPM particle-file-writer cases.

This script is intentionally stdlib-only.  It prepares isolated run directories,
invokes particleFileWriter.py, optionally submits the generated GEOS batch job,
and writes a Markdown/JSON report that can be collected from Dane after jobs run.
"""
from __future__ import annotations

import argparse
import ast
import datetime as _dt
import fnmatch
import getpass
import json
import os
import re
import shutil
import subprocess
import sys
import textwrap
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

PFW_INPUT_PREFIX = "pfw_input_"
STATUS_FILE = "pfw_suite_status.json"
DEFAULT_CORE_FILES = (
    "particleFileWriter.py",
    "pfw_check.py",
    "pfw_geometryObjects.py",
    "pfw_materials.py",
    "pfw_analysis.py",
    "pfw_visit_render.py",
)
DEFAULT_EXCLUDE_DIRS = {
    "__pycache__",
    ".git",
    ".pytest_cache",
    "_suite_runs",
    "suite_runs",
    "suite_results",
    "visit_frames",
}
TEMPLATE_TOKEN_RE = re.compile(r"(XXXX|YYYY|ZZZZ|<[^>]+>)")
JOB_ID_RE = re.compile(r"Submitted job with ID\s*=\s*(\S+)|Submitted batch job\s+(\S+)")


@dataclass(frozen=True)
class Case:
    case_id: str
    input_path: Path
    suite_root: Path
    description: str
    is_template: bool = False

    @property
    def input_name(self) -> str:
        return self.input_path.name

    @property
    def input_module(self) -> str:
        return self.input_path.stem

    @property
    def short_name(self) -> str:
        stem = self.input_path.stem
        if stem.startswith(PFW_INPUT_PREFIX):
            return stem[len(PFW_INPUT_PREFIX):]
        return stem

    @property
    def source_dir(self) -> Path:
        return self.input_path.parent


@dataclass
class CaseResult:
    case_id: str
    source: str
    run_dir: str
    status: str = "not_started"
    stage: str = "not_started"
    returncode: int | None = None
    job_id: str | None = None
    submit_requested: bool = False
    generated_xml: list[str] = field(default_factory=list)
    generated_particle_files: list[str] = field(default_factory=list)
    reaction_history: bool = False
    box_average_history: bool = False
    visit_frames: list[str] = field(default_factory=list)
    message: str = ""
    started: str | None = None
    finished: str | None = None


def now_iso() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def pfw_root() -> Path:
    return Path(__file__).resolve().parent


def relpath(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def normalize_path(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()


def extract_description(path: Path, max_lines: int = 16) -> str:
    """Return the leading comment block from a pfw_input file."""
    lines: list[str] = []
    try:
        raw_lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return ""

    for line in raw_lines:
        stripped = line.strip()
        if stripped.startswith("# -*-"):
            continue
        if not stripped:
            if lines:
                lines.append("")
            continue
        if stripped.startswith("#"):
            text = stripped.lstrip("#").strip()
            if text and not set(text) <= {"-", "=", "_"}:
                lines.append(text)
            continue
        break

    compact: list[str] = []
    previous_blank = False
    for line in lines:
        blank = not line.strip()
        if blank and previous_blank:
            continue
        compact.append(line)
        previous_blank = blank

    return " ".join(x for x in compact[:max_lines] if x).strip()


def discover_cases(
    suite_root: Path,
    recursive: bool = True,
    include_templates: bool = False,
    only: Sequence[str] | None = None,
    exclude: Sequence[str] | None = None,
) -> list[Case]:
    suite_root = normalize_path(suite_root)
    globber = suite_root.rglob if recursive else suite_root.glob
    only = list(only or [])
    exclude = list(exclude or [])

    cases: list[Case] = []
    for path in sorted(globber(f"{PFW_INPUT_PREFIX}*.py")):
        if any(part in DEFAULT_EXCLUDE_DIRS for part in path.parts):
            continue
        if path.name.endswith("_suite.py"):
            continue
        rel = path.relative_to(suite_root)
        stem = path.stem[len(PFW_INPUT_PREFIX):]
        parent = rel.parent
        case_id = stem if str(parent) == "." else (parent / stem).as_posix()
        is_template = bool(TEMPLATE_TOKEN_RE.search(path.name) or TEMPLATE_TOKEN_RE.search(path.read_text(errors="replace")[:512]))

        if is_template and not include_templates:
            continue
        if only and not any(fnmatch.fnmatch(case_id, pat) or fnmatch.fnmatch(path.as_posix(), pat) for pat in only):
            continue
        if exclude and any(fnmatch.fnmatch(case_id, pat) or fnmatch.fnmatch(path.as_posix(), pat) for pat in exclude):
            continue

        cases.append(Case(case_id=case_id, input_path=path, suite_root=suite_root,
                          description=extract_description(path), is_template=is_template))
    return cases


def compile_preflight(case: Case) -> dict:
    result: dict = {
        "case_id": case.case_id,
        "source": relpath(case.input_path, case.suite_root),
        "ok": True,
        "errors": [],
        "warnings": [],
        "description": case.description,
        "is_template": case.is_template,
    }
    text = case.input_path.read_text(errors="replace")
    try:
        tree = ast.parse(text, filename=str(case.input_path))
    except SyntaxError as exc:
        result["ok"] = False
        result["errors"].append(f"SyntaxError line {exc.lineno}: {exc.msg}")
        return result

    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            for keyword in node.keywords:
                if keyword.arg == "v":
                    line = getattr(node, "lineno", "?")
                    result["warnings"].append(f"Deprecated geometry keyword v= in call starting at line {line}; use vel=")
    return result


def print_preflight(results: list[dict]) -> int:
    failures = [r for r in results if not r["ok"]]
    warnings = [(r, w) for r in results for w in r["warnings"]]
    print(f"Preflight cases: {len(results)}")
    print(f"  errors:   {len(failures)}")
    print(f"  warnings: {len(warnings)}")
    if failures:
        print("\nErrors:")
        for r in failures:
            for e in r["errors"]:
                print(f"  - {r['case_id']}: {e}")
    if warnings:
        print("\nWarnings:")
        for r, w in warnings:
            print(f"  - {r['case_id']}: {w}")
    return 1 if failures else 0


def copy_case_files(case: Case, run_dir: Path) -> None:
    """Copy the source problem directory into run_dir, skipping transient output."""
    run_dir.mkdir(parents=True, exist_ok=True)
    for src in case.source_dir.iterdir():
        if src.name in DEFAULT_EXCLUDE_DIRS:
            continue
        dst = run_dir / src.name
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*DEFAULT_EXCLUDE_DIRS))
        else:
            shutil.copy2(src, dst)


def copy_core_files(root: Path, run_dir: Path) -> None:
    for name in DEFAULT_CORE_FILES:
        src = root / name
        if src.exists():
            shutil.copy2(src, run_dir / name)


def write_user_defs(run_dir: Path, root: Path, geos_path: str, run_root: Path, bank: str | None) -> Path:
    username = getpass.getuser()
    user_defs = run_dir / f"userDefs_{username}.py"
    default_bank = bank or os.environ.get("GEOS_BANK", "")
    user_defs.write_text(textwrap.dedent(f"""
        # Auto-generated by pfw_suite.py.  Do not commit machine-local paths.
        import platform

        pfwPath = r"{root.as_posix()}"
        geosPath = r"{geos_path}"
        testRunDirectory = r"{run_root.as_posix()}"
        defaultRunDirectory = r"{run_root.as_posix()}"
        defaultBank = r"{default_bank}"
        lassen = 'lassen' in platform.node()
        dane = 'dane' in platform.node()
        tuolumne = 'tuolumne' in platform.node()
        rzhound = 'rzhound' in platform.node()
        tioga = 'tioga' in platform.node()
        """).lstrip())
    return user_defs


def append_suite_overrides(
    original_input: Path,
    run_input: Path,
    submit: bool,
    bank: str | None,
    partition: str | None,
    wall_time: str | None,
    output_type: str | None,
    run_debug: bool | None,
    auto_restart: bool | None,
) -> None:
    text = original_input.read_text(errors="replace")
    overrides: list[str] = [
        "",
        "# ---- pfw_suite.py overrides -------------------------------------------------",
        f"pfw['mSubmitJobs'] = {bool(submit)!r}",
    ]
    if bank:
        overrides.append(f"pfw['mBank'] = {bank!r}")
    if partition:
        overrides.append(f"pfw['mPartition'] = {partition!r}")
    if wall_time:
        overrides.append(f"pfw['mWallTime'] = {wall_time!r}")
    if output_type:
        overrides.append(f"pfw['outputType'] = {output_type!r}")
    if run_debug is not None:
        overrides.append(f"pfw['runDebug'] = {bool(run_debug)!r}")
    if auto_restart is not None:
        overrides.append(f"pfw['autoRestart'] = {bool(auto_restart)!r}")
    overrides.append("# -----------------------------------------------------------------------------")
    run_input.write_text(text.rstrip() + "\n" + "\n".join(overrides) + "\n")


def run_dir_for_case(run_root: Path, suite_name: str, case: Case) -> Path:
    return run_root / suite_name / case.case_id


def safe_remove(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)


def parse_job_id(stdout: str) -> str | None:
    for match in JOB_ID_RE.finditer(stdout):
        return next(group for group in match.groups() if group)
    return None


def collect_generated_files(run_dir: Path) -> tuple[list[str], list[str]]:
    xmls = sorted(p.name for p in run_dir.glob("mpm_*.xml"))
    particles = sorted(p.name for p in run_dir.glob("mpmParticleFile_*"))
    return xmls, particles


def write_status(run_dir: Path, result: CaseResult) -> None:
    (run_dir / STATUS_FILE).write_text(json.dumps(asdict(result), indent=2, sort_keys=True) + "\n")


def load_status(run_dir: Path, case: Case | None = None) -> CaseResult | None:
    path = run_dir / STATUS_FILE
    if not path.exists():
        if case is None:
            return None
        return CaseResult(case_id=case.case_id, source=str(case.input_path), run_dir=str(run_dir))
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    return CaseResult(**data)


def prepare_case(case: Case, args: argparse.Namespace, root: Path) -> CaseResult:
    run_dir = run_dir_for_case(args.run_root, args.suite_name, case)
    result = CaseResult(case_id=case.case_id, source=str(case.input_path), run_dir=str(run_dir),
                        status="prepared", stage="prepare", submit_requested=args.submit)
    if args.clean:
        safe_remove(run_dir)
    if run_dir.exists() and any(run_dir.iterdir()) and not args.force and not args.clean:
        result.status = "skipped"
        result.message = "Run directory exists; use --force or --clean to overwrite prepared inputs."
        return result

    run_dir.mkdir(parents=True, exist_ok=True)
    copy_case_files(case, run_dir)
    copy_core_files(root, run_dir)
    write_user_defs(run_dir, root, args.geos_path, args.run_root, args.bank)
    append_suite_overrides(
        case.input_path,
        run_dir / case.input_name,
        submit=args.submit,
        bank=args.bank,
        partition=args.partition,
        wall_time=args.wall_time,
        output_type=args.output_type,
        run_debug=args.run_debug,
        auto_restart=args.auto_restart,
    )
    write_status(run_dir, result)
    return result


def invoke_pfw(case: Case, args: argparse.Namespace, result: CaseResult) -> CaseResult:
    run_dir = Path(result.run_dir)
    result.stage = "particle_file_writer"
    result.status = "running"
    result.started = now_iso()
    write_status(run_dir, result)

    python_exe = args.python or sys.executable or "python3"
    cmd = [python_exe, "particleFileWriter.py", case.input_module]
    if args.pfw_ranks and args.pfw_ranks > 1:
        cmd = ["srun", "-n", str(args.pfw_ranks)] + cmd

    stdout_path = run_dir / "pfw_generate.out"
    stderr_path = run_dir / "pfw_generate.err"
    if args.dry_run:
        result.status = "dry_run"
        result.stage = "prepared"
        result.message = "Prepared run directory; did not invoke particleFileWriter.py."
        result.finished = now_iso()
        write_status(run_dir, result)
        return result

    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        proc = subprocess.run(cmd, cwd=run_dir, stdout=stdout, stderr=stderr, text=True)

    result.returncode = proc.returncode
    result.finished = now_iso()
    stdout_text = stdout_path.read_text(errors="replace") if stdout_path.exists() else ""
    result.job_id = parse_job_id(stdout_text)
    result.generated_xml, result.generated_particle_files = collect_generated_files(run_dir)
    if proc.returncode == 0:
        result.status = "submitted" if result.job_id else "generated"
        result.stage = "submitted" if result.job_id else "generated"
    else:
        result.status = "failed"
        result.stage = "particle_file_writer_failed"
        tail = "\n".join((stderr_path.read_text(errors="replace") + stdout_text).splitlines()[-20:])
        result.message = tail
    write_status(run_dir, result)
    return result


def infer_run_state(run_dir: Path, result: CaseResult | None = None) -> CaseResult:
    if result is None:
        result = CaseResult(case_id=run_dir.name, source="", run_dir=str(run_dir))

    result.generated_xml, result.generated_particle_files = collect_generated_files(run_dir)
    result.reaction_history = (run_dir / "reactionHistory.csv").exists()
    result.box_average_history = (run_dir / "boxAverageHistory.csv").exists()
    frame_dir = run_dir / "visit_frames"
    if frame_dir.exists():
        result.visit_frames = sorted(relpath(p, run_dir) for p in frame_dir.glob("*.png"))

    slurm_logs = sorted(run_dir.glob("slurm-*.out"), key=lambda p: p.stat().st_mtime)
    if slurm_logs:
        log_text = slurm_logs[-1].read_text(errors="replace")
        if "Job complete" in log_text or "srun command has completed" in log_text or "flux run" in log_text and "completed" in log_text:
            result.status = "complete"
            result.stage = "geos_complete"
        elif "TIME LIMIT" in log_text or "NODE FAILURE" in log_text or "Job exited early" in log_text:
            result.status = "needs_restart_or_failed"
            result.stage = "geos_interrupted"
        elif "error" in log_text.lower() or "segmentation" in log_text.lower():
            result.status = "failed"
            result.stage = "geos_failed"
        elif result.status in {"submitted", "not_started"}:
            result.status = "submitted_or_running"
            result.stage = "geos_submitted"
    elif result.generated_xml and result.generated_particle_files and result.status == "not_started":
        result.status = "generated"
        result.stage = "generated"

    return result


def write_markdown_report(
    cases: list[Case],
    results: list[CaseResult],
    output: Path,
    suite_name: str,
    run_root: Path,
    geos_path: str | None = None,
) -> None:
    by_id = {r.case_id: r for r in results}
    counts: dict[str, int] = {}
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1

    lines = []
    lines.append(f"# GEOS-MPM {suite_name} suite report")
    lines.append("")
    lines.append(f"Generated: {now_iso()}")
    lines.append(f"Run root: `{run_root}`")
    if geos_path:
        lines.append(f"GEOS executable: `{geos_path}`")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"Cases discovered: **{len(cases)}**")
    for status in sorted(counts):
        lines.append(f"- {status}: {counts[status]}")
    lines.append("")
    lines.append("## Cases")
    lines.append("")
    lines.append("| Case | Status | Stage | Job ID | XML | Particle file | Reactions | Box averages | Description |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for case in cases:
        r = by_id.get(case.case_id) or CaseResult(case_id=case.case_id, source=str(case.input_path),
                                                  run_dir=str(run_dir_for_case(run_root, suite_name, case)))
        xml = ", ".join(r.generated_xml)
        pf = ", ".join(r.generated_particle_files)
        desc = case.description.replace("|", "\\|") if case.description else ""
        lines.append(
            f"| `{case.case_id}` | {r.status} | {r.stage} | {r.job_id or ''} | {xml} | {pf} | "
            f"{'yes' if r.reaction_history else ''} | {'yes' if r.box_average_history else ''} | {desc} |"
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("`generated` means pfw completed and wrote the GEOS XML/particle file but did not submit a batch job. `submitted` means pfw reported a scheduler job id. `complete` is inferred from scheduler output in the run directory.")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")


def write_json_report(results: list[CaseResult], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps([asdict(r) for r in results], indent=2, sort_keys=True) + "\n")


def run_cases(args: argparse.Namespace) -> int:
    root = pfw_root()
    cases = discover_cases(args.suite_root, recursive=args.recursive,
                           include_templates=args.include_templates, only=args.only, exclude=args.exclude)
    if not cases:
        print("No cases matched.", file=sys.stderr)
        return 2

    preflight = [compile_preflight(c) for c in cases]
    if args.preflight:
        code = print_preflight(preflight)
        if code and not args.keep_going:
            return code

    if args.list:
        for case in cases:
            template = " [template]" if case.is_template else ""
            print(f"{case.case_id}{template}\t{relpath(case.input_path, args.suite_root)}")
        return 0

    results: list[CaseResult] = []
    for case, pf in zip(cases, preflight):
        if not pf["ok"]:
            run_dir = run_dir_for_case(args.run_root, args.suite_name, case)
            r = CaseResult(case_id=case.case_id, source=str(case.input_path), run_dir=str(run_dir),
                           status="preflight_failed", stage="preflight",
                           message="; ".join(pf["errors"]))
            results.append(r)
            print(f"[preflight failed] {case.case_id}: {r.message}")
            if not args.keep_going:
                break
            continue

        print(f"[prepare] {case.case_id}")
        result = prepare_case(case, args, root)
        if result.status == "skipped":
            print(f"[skip] {case.case_id}: {result.message}")
            results.append(result)
            continue
        print(f"[pfw] {case.case_id}")
        result = invoke_pfw(case, args, result)
        print(f"[{result.status}] {case.case_id}")
        results.append(result)
        if result.status == "failed" and not args.keep_going:
            break

    report_path = args.report or (args.run_root / args.suite_name / "suite_report.md")
    json_path = args.json_report or (args.run_root / args.suite_name / "suite_report.json")
    status_results = [infer_run_state(Path(r.run_dir), r) for r in results]
    write_markdown_report(cases, status_results, report_path, args.suite_name, args.run_root, args.geos_path)
    write_json_report(status_results, json_path)
    print(f"Report: {report_path}")
    print(f"JSON:   {json_path}")
    return 0 if all(r.status not in {"failed", "preflight_failed"} for r in status_results) else 1


def status_cases(args: argparse.Namespace) -> int:
    cases = discover_cases(args.suite_root, recursive=args.recursive,
                           include_templates=args.include_templates, only=args.only, exclude=args.exclude)
    results: list[CaseResult] = []
    for case in cases:
        run_dir = run_dir_for_case(args.run_root, args.suite_name, case)
        result = load_status(run_dir, case)
        if result is None:
            result = CaseResult(case_id=case.case_id, source=str(case.input_path), run_dir=str(run_dir))
        results.append(infer_run_state(run_dir, result))
    report_path = args.report or (args.run_root / args.suite_name / "suite_report.md")
    json_path = args.json_report or (args.run_root / args.suite_name / "suite_report.json")
    write_markdown_report(cases, results, report_path, args.suite_name, args.run_root, args.geos_path)
    write_json_report(results, json_path)
    for r in results:
        print(f"{r.status:24s} {r.case_id}")
    print(f"Report: {report_path}")
    return 0


def emit_run_script(directory: Path, recursive: bool) -> None:
    script = directory / "runProblem"
    recursive_flag = "" if recursive else " --no-recursive"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        set -euo pipefail
        SCRIPT_DIR=$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd)
        PFW_ROOT="$SCRIPT_DIR"
        while [ ! -f "$PFW_ROOT/pfw_suite.py" ] && [ "$PFW_ROOT" != "/" ]; do
          PFW_ROOT=$(dirname "$PFW_ROOT")
        done
        if [ ! -f "$PFW_ROOT/pfw_suite.py" ]; then
          echo "Could not locate pfw_suite.py" >&2
          exit 2
        fi
        action="${{1:-run}}"
        case "$action" in
          run|list|preflight|status|report)
            shift || true
            exec python3 "$PFW_ROOT/pfw_suite.py" "$action" --suite-root "$SCRIPT_DIR" --suite-name "$(basename "$SCRIPT_DIR")"{recursive_flag} "$@"
            ;;
          *)
            exec python3 "$PFW_ROOT/pfw_suite.py" run --suite-root "$SCRIPT_DIR" --suite-name "$(basename "$SCRIPT_DIR")"{recursive_flag} "$@"
            ;;
        esac
        """))
    script.chmod(0o755)


def emit_run_scripts(args: argparse.Namespace) -> int:
    suite_root = normalize_path(args.suite_root)
    dirs = sorted({p.parent for p in suite_root.rglob(f"{PFW_INPUT_PREFIX}*.py")})
    count = 0
    for directory in dirs:
        if any(part in DEFAULT_EXCLUDE_DIRS for part in directory.parts):
            continue
        emit_run_script(directory, recursive=args.problem_recursive)
        count += 1
    print(f"Wrote {count} runProblem scripts under {suite_root}")
    return 0


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--suite-root", type=normalize_path, default=Path.cwd(), help="Directory containing pfw_input_*.py files")
    parser.add_argument("--suite-name", default=None, help="Name used below --run-root; defaults to suite-root basename")
    parser.add_argument("--run-root", type=normalize_path, default=None, help="Root directory for prepared/running jobs")
    parser.add_argument("--recursive", dest="recursive", action="store_true", default=True, help="Discover cases recursively")
    parser.add_argument("--no-recursive", dest="recursive", action="store_false", help="Only discover cases directly in --suite-root")
    parser.add_argument("--include-templates", action="store_true", help="Include template inputs containing XXXX/YYYY/ZZZZ tokens")
    parser.add_argument("--only", action="append", default=[], help="fnmatch pattern for case ids to include; may be repeated")
    parser.add_argument("--exclude", action="append", default=[], help="fnmatch pattern for case ids to exclude; may be repeated")
    parser.add_argument("--geos-path", default=os.environ.get("GEOS_EXECUTABLE", "geosx"), help="Path to GEOS executable written into userDefs")
    parser.add_argument("--bank", default=os.environ.get("GEOS_BANK"), help="LC bank/account for generated scripts")
    parser.add_argument("--partition", default=None, help="Scheduler partition/queue override, e.g. pdebug or pbatch")
    parser.add_argument("--wall-time", default=None, help="Walltime override, HH:MM:SS")
    parser.add_argument("--output-type", choices=("silo", "vtk"), default=os.environ.get("PFW_OUTPUT_TYPE", "silo"), help="Override pfw['outputType']; defaults to silo for the Dane MPM minimal-TPL build")
    parser.add_argument("--report", type=normalize_path, default=None, help="Markdown report path")
    parser.add_argument("--json-report", type=normalize_path, default=None, help="JSON report path")


def finalize_common_args(args: argparse.Namespace) -> argparse.Namespace:
    args.suite_root = normalize_path(args.suite_root)
    if args.suite_name is None:
        args.suite_name = args.suite_root.name
    if args.run_root is None:
        args.run_root = args.suite_root / "_suite_runs"
    else:
        args.run_root = normalize_path(args.run_root)
    return args


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run and report GEOS-MPM pfw suites")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Prepare run directories and invoke particleFileWriter.py")
    add_common_args(run)
    run.add_argument("--submit", action="store_true", help="Allow pfw to submit generated GEOS batch jobs")
    run.add_argument("--no-submit", dest="submit", action="store_false", help="Generate only; do not submit jobs")
    run.set_defaults(submit=False)
    run.add_argument("--clean", action="store_true", help="Delete existing run directories before preparing")
    run.add_argument("--force", action="store_true", help="Overwrite prepared files in existing run directories")
    run.add_argument("--dry-run", action="store_true", help="Prepare run directories but do not invoke pfw")
    run.add_argument("--preflight", action="store_true", default=True, help="Compile-check inputs before running")
    run.add_argument("--no-preflight", dest="preflight", action="store_false")
    run.add_argument("--keep-going", action="store_true", help="Continue after a preflight/pfw failure")
    run.add_argument("--pfw-ranks", type=int, default=1, help="Use srun -n N for parallel particle generation")
    run.add_argument("--python", default=None, help="Python executable for particleFileWriter.py")
    debug_group = run.add_mutually_exclusive_group()
    debug_group.add_argument("--debug", dest="run_debug", action="store_true", help="Set pfw['runDebug']=True")
    debug_group.add_argument("--release-queue", dest="run_debug", action="store_false", help="Set pfw['runDebug']=False")
    run.set_defaults(run_debug=None)
    restart_group = run.add_mutually_exclusive_group()
    restart_group.add_argument("--auto-restart", dest="auto_restart", action="store_true", help="Set pfw['autoRestart']=True")
    restart_group.add_argument("--no-auto-restart", dest="auto_restart", action="store_false", help="Set pfw['autoRestart']=False")
    run.set_defaults(auto_restart=None)
    run.add_argument("--list", action="store_true", help="List matching cases and exit")

    for name, help_text in (("list", "List cases"), ("preflight", "Compile-check inputs"), ("status", "Scan run directories and write report"), ("report", "Alias for status")):
        p = sub.add_parser(name, help=help_text)
        add_common_args(p)

    emit = sub.add_parser("emit-run-scripts", help="Create runProblem wrappers in problem directories")
    add_common_args(emit)
    emit.add_argument("--problem-recursive", action="store_true", help="Generated runProblem wrappers also discover nested child cases")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args = finalize_common_args(args)

    if args.command == "run":
        return run_cases(args)
    if args.command == "list":
        cases = discover_cases(args.suite_root, recursive=args.recursive,
                               include_templates=args.include_templates, only=args.only, exclude=args.exclude)
        for case in cases:
            print(f"{case.case_id}\t{relpath(case.input_path, args.suite_root)}")
        return 0
    if args.command == "preflight":
        cases = discover_cases(args.suite_root, recursive=args.recursive,
                               include_templates=args.include_templates, only=args.only, exclude=args.exclude)
        return print_preflight([compile_preflight(c) for c in cases])
    if args.command in {"status", "report"}:
        return status_cases(args)
    if args.command == "emit-run-scripts":
        return emit_run_scripts(args)
    parser.error(f"Unknown command {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
