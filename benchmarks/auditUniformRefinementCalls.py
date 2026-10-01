#!/usr/bin/env python3
"""Verify runtime call counts in the actual positive uniform-refinement interval."""

import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--expect-probe", action="store_true", help="require injected forbidden calls as an instrumentation self-test")
    parser.add_argument("--minimum-calls", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.minimum_calls < 1 or args.timeout <= 0:
        parser.error("a command and positive minimum call count are required")
    env = os.environ.copy()
    previous = env.get("LD_PRELOAD", "")
    env["LD_PRELOAD"] = str(args.library.resolve(strict=True)) + (":" + previous if previous else "")
    if args.expect_probe:
        env["GEOS_REFINEMENT_AUDIT_PROBE"] = "1"
    else:
        env.pop("GEOS_REFINEMENT_AUDIT_PROBE", None)
    process = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, start_new_session=True)
    try:
        output, _ = process.communicate(timeout=args.timeout)
    except (subprocess.TimeoutExpired, KeyboardInterrupt):
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            output, _ = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
        if args.log:
            args.log.write_text(output)
        raise
    if args.log:
        args.log.write_text(output)
    records = []
    for match in re.finditer(r"uniform_refinement_audit (rank=\d+ levels=\d+ mpi=\d+ gather=\d+ scatter=\d+ redistribute=\d+ append=\d+ clean=\d+ repartition=\d+)", output):
        records.append({key: int(value) for key, value in (pair.split("=") for pair in match.group(1).split())})
    error = None
    if process.returncode:
        error = f"command failed with status {process.returncode}"
    elif len(records) < args.minimum_calls:
        error = f"only {len(records)} observed refinement calls; require {args.minimum_calls}"
    elif args.expect_probe:
        if any(record["append"] < 1 or record["clean"] < 1 or
               (record["mpi"] and (record["gather"] < 1 or record["repartition"] < 1)) for record in records):
            error = "injected forbidden calls were not observed"
    elif any(record[key] for record in records for key in ("gather", "scatter", "redistribute", "append", "clean", "repartition")):
        error = "forbidden work occurred inside positive refinement"
    print(json.dumps({"command": command, "probe": args.expect_probe, "calls": records, "error": error}, indent=2))
    if error:
        if not args.log:
            print(output, file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
