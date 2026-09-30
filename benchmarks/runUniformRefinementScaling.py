#!/usr/bin/env python3
"""Run component strong/weak scaling inside an existing MPI/Slurm allocation."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import shlex
import statistics
import subprocess
import sys
import time


def positive(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def provenance(binary, repo, discovery, sharing):
    def digest(path):
        hasher = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                hasher.update(chunk)
        return hasher.hexdigest()

    def git(*args):
        result = subprocess.run(
            ["git", "-C", str(repo), *args], capture_output=True, text=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else None

    files = [
        *repo.glob("src/coreComponents/mesh/generators/VTKRefinement*.[ch]pp"),
        *repo.glob("src/coreComponents/common/MpiChunkedCommunication*.[ch]pp"),
        *repo.glob("src/coreComponents/mesh/generators/VTKMeshScattering*.[ch]pp"),
        repo / "src/coreComponents/mesh/benchmarks/benchmarkUniformRefinement.cpp",
        Path(__file__).resolve(),
    ]
    cache = binary.parent.parent / "CMakeCache.txt"
    keys = (
        "CMAKE_BUILD_TYPE:", "CMAKE_CXX_COMPILER:", "CMAKE_CXX_FLAGS:",
        "GEOS_TPL_DIR:", "ENABLE_MPI:", "ENABLE_OPENMP:", "ENABLE_CUDA:",
        "ENABLE_HIP:", "VTK_DIR:", "MPI_CXX_COMPILER:",
    )
    return {
        "component_only": True,
        "discovery": "boundary geometry plus all volume IDs" if discovery == "boundary" else "all coarse entities (diagnostic fallback)",
        "incidence": "shared interfaces only" if sharing == "interfaces" else "all local entities on the current level",
        "host": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "binary": str(binary),
        "binary_sha256": digest(binary),
        "component_library_sha256": {
            str(path): digest(path)
            for pattern in ("libvtkRefinement*.so", "libcommon.so", "libmesh.so", "liblvarray.so")
            for path in (binary.parent.parent / "lib").glob(pattern) if path.is_file()
        },
        "git_head": git("rev-parse", "HEAD"),
        "git_status": git("status", "--short"),
        "source_sha256": {
            str(path.relative_to(repo)): digest(path)
            for path in files if path.is_file()
        },
        "build_cache": [
            line for line in cache.read_text().splitlines() if line.startswith(keys)
        ] if cache.is_file() else [],
        "environment": {
            key: os.environ[key] for key in (
                "SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_CPUS_PER_TASK",
                "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "LOADEDMODULES",
                "VTK_SMP_MAX_THREADS",
            ) if key in os.environ
        },
    }


def summarize(directory):
    groups = {}
    for path in sorted(directory.glob("*.csv")):
        if path.name == "summary.csv":
            continue
        metadata_path = path.with_suffix(".json")
        if not metadata_path.is_file():
            continue
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("returncode") != 0:
            continue
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        totals = [row for row in rows if row["phase"] == "levelTotal"]
        if not totals:
            continue
        key = tuple(metadata[name] for name in ("mode", "kind", "nodes", "ranks", "levels"))
        # Sum per-level critical-rank times; reporting and phase barriers are excluded.
        groups.setdefault(key, []).append({
            "seconds": sum(float(row["seconds_max"]) for row in totals),
            "cells": int(totals[-1]["cells_sum"]),
            "imbalance": int(totals[-1]["cells_max"]) * metadata["ranks"] / int(totals[-1]["cells_sum"]),
            "peak_rss_kib": max(int(row["peak_rss_kib_max"]) for row in rows),
            "payload_bytes": sum(int(row["payload_bytes_sum"]) for row in totals),
        })
    records = []
    for key, repeats in sorted(groups.items()):
        mode, kind, nodes, ranks, levels = key
        baseline_key = min(
            (candidate for candidate in groups if candidate[0] == mode and candidate[1] == kind and candidate[4] == levels),
            key=lambda candidate: candidate[3],
        )
        baseline = statistics.median(item["seconds"] for item in groups[baseline_key])
        seconds = statistics.median(item["seconds"] for item in repeats)
        efficiency = baseline / seconds
        if mode == "strong":
            efficiency *= baseline_key[3] / ranks
        records.append({
            "mode": mode, "kind": kind, "nodes": nodes, "ranks": ranks,
            "levels": levels, "repeats": len(repeats),
            "median_refinement_seconds": seconds,
            "min_refinement_seconds": min(item["seconds"] for item in repeats),
            "max_refinement_seconds": max(item["seconds"] for item in repeats),
            "baseline_ranks": baseline_key[3], "efficiency": efficiency,
            "fine_cells": repeats[0]["cells"],
            "cells_per_second": repeats[0]["cells"] / seconds,
            "cell_imbalance": max(item["imbalance"] for item in repeats),
            "peak_rank_rss_kib": max(item["peak_rss_kib"] for item in repeats),
            "refinement_payload_bytes": statistics.median(item["payload_bytes"] for item in repeats),
        })
    if records:
        with (directory / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path, help="Release/RelWithDebInfo benchmarkUniformRefinement")
    parser.add_argument("output", type=Path, help="new results directory")
    parser.add_argument("--launcher", default="mpiexec", help='e.g. "srun --cpu-bind=cores"')
    parser.add_argument("--nodes", type=positive, nargs="+", default=[1, 2, 4])
    parser.add_argument("--ranks-per-node", type=positive, required=True)
    parser.add_argument("--strong-cells", type=positive, nargs=3, default=[32, 32, 32])
    parser.add_argument("--weak-cells", type=positive, nargs=3, default=[8, 8, 8])
    parser.add_argument("--levels", type=positive, nargs="+", default=[1, 2])
    parser.add_argument("--kinds", choices=["hex", "tet", "pyramid"], nargs="+", default=["hex", "tet", "pyramid"])
    parser.add_argument("--modes", choices=["strong", "weak"], nargs="+", default=["strong", "weak"])
    parser.add_argument("--repeats", type=positive, default=3)
    parser.add_argument("--point-components", type=int, default=0)
    parser.add_argument("--discovery", choices=["boundary", "all"], default="boundary")
    parser.add_argument("--sharing", choices=["interfaces", "all"], default="interfaces")
    parser.add_argument("--input", type=Path, help="VTU for strong scaling only; uses Cartesian coarse scatter")
    parser.add_argument("--timeout", type=positive, default=1800, help="seconds per MPI invocation")
    parser.add_argument("--dry-run", action="store_true", help="write commands/metadata without launching MPI")
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        parser.error("binary must be an executable file")
    if args.point_components < 0:
        parser.error("--point-components must be nonnegative")
    if args.input and (args.modes != ["strong"] or args.point_components):
        parser.error("--input requires --modes strong and --point-components 0")
    launcher = shlex.split(args.launcher)
    if not launcher:
        parser.error("empty launcher")
    slurm = Path(launcher[0]).name == "srun"
    if slurm and not args.dry_run and "SLURM_JOB_ID" not in os.environ:
        parser.error("run srun cases inside an existing Slurm allocation")
    args.output.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parent.parent
    base = provenance(binary, repo, args.discovery, args.sharing)
    (args.output / "provenance.json").write_text(json.dumps(base, indent=2) + "\n")
    cases = []
    for mode in args.modes:
        for kind in (["input"] if args.input else args.kinds):
            for nodes in sorted(set(args.nodes)):
                ranks = nodes * args.ranks_per_node
                for levels in sorted(set(args.levels)):
                    for repeat in range(args.repeats):
                        label = f"{mode}-{kind}-n{nodes}-p{ranks}-l{levels}-r{repeat}"
                        command = launcher + (
                            [f"--nodes={nodes}", f"--ntasks={ranks}", f"--ntasks-per-node={args.ranks_per_node}"]
                            if slurm else ["-n", str(ranks)]
                        ) + [str(binary), "--levels", str(levels), "--label", label,
                             "--discovery", args.discovery, "--sharing", args.sharing]
                        if args.input:
                            command += ["--input", str(args.input.resolve())]
                        else:
                            cells = args.weak_cells if mode == "weak" else args.strong_cells
                            command += ["--kind", kind, "--cells", *map(str, cells), "--point-components", str(args.point_components)]
                        if mode == "weak":
                            command.append("--weak")
                        cases.append({
                            "label": label, "mode": mode, "kind": kind, "nodes": nodes,
                            "ranks": ranks, "levels": levels, "repeat": repeat, "command": command,
                        })
    (args.output / "commands.json").write_text(json.dumps(cases, indent=2) + "\n")
    for case in cases:
        print(shlex.join(case["command"]), flush=True)
        if args.dry_run:
            continue
        case["start_unix_seconds"] = time.time()
        output = args.output / case["label"]
        try:
            with output.with_suffix(".csv").open("w") as stdout, output.with_suffix(".stderr").open("w") as stderr:
                process = subprocess.Popen(case["command"], stdout=stdout, stderr=stderr, start_new_session=True)
                try:
                    case["returncode"] = process.wait(timeout=args.timeout)
                except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                    case["returncode"] = "timeout" if isinstance(error, subprocess.TimeoutExpired) else "interrupted"
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
        except OSError as error:
            case["returncode"] = str(error)
        case["wall_seconds"] = time.time() - case["start_unix_seconds"]
        output.with_suffix(".json").write_text(json.dumps(case, indent=2) + "\n")
        if case["returncode"] != 0:
            summarize(args.output)
            raise SystemExit(f"Benchmark failed ({case['returncode']}): {output}.stderr")
    summarize(args.output)


if __name__ == "__main__":
    main()
