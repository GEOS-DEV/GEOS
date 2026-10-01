#!/usr/bin/env python3
"""Measure the full GEOS refinement/import/ghost/DOF/solve path in an allocation.

Generated meshes are connected hex grids with a pressure-gradient FVM solve.
--mesh supplies an existing VTU (including mixed cells/prisms) or VTM; that case
uses a uniform nonzero source with closed boundaries and supports strong scaling.
--face-blocks adds fracture geometry/relations/ghost initialization and coupled
matrix/fracture single-phase flow with a fixed aperture to a VTM run.
Component measurements are provided separately by runUniformRefinementScaling.py.
"""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import socket
import statistics
import subprocess
import sys
from xml.etree import ElementTree as ET

from runUniformRefinementScaling import positive, provenance
from verifyUniformRefinementInitialization import make_inputs


def rank_wrapper(arguments):
    """Run on every rank; record actual placement and GNU time's peak RSS."""
    directory, executable, *command = arguments
    rank = next((os.environ[name] for name in (
        "PMI_RANK", "PMIX_RANK", "OMPI_COMM_WORLD_RANK", "SLURM_PROCID")
        if name in os.environ), "0")
    path = Path(directory) / f"rank-{int(rank):07}"
    metadata = {"rank": int(rank), "host": socket.gethostname(),
                "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None}
    path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    os.execv("/usr/bin/time", ["/usr/bin/time", "-f", "%e,%M,%x", "-o",
                              str(path.with_suffix(".time")), executable, *command])


def process_grid(ranks):
    triples = [(a, b, ranks // (a * b)) for a in range(1, math.isqrt(ranks) + 1)
               if ranks % a == 0 for b in range(a, math.isqrt(ranks // a) + 1)
               if (ranks // a) % b == 0]
    return tuple(reversed(min(triples, key=lambda t: t[0] * t[1] + t[1] * t[2] + t[2] * t[0])))


def make_case_inputs(directory, cells, domain, mesh, partition_method, partition_refinement,
                     main_block="main", face_blocks=()):
    directory.mkdir()
    make_inputs(directory, cells, domain)
    paths = {}
    for level in (0, 1, 2):
        path = directory / f"fvm-L{level}.xml"
        tree = ET.parse(path)
        root = tree.getroot()
        root.find("Events").set("maxTime", "1")
        events = root.find("Events")
        for event in list(events):
            if event.get("target", "").startswith("/Outputs/"):
                events.remove(event)
        root.remove(root.find("Outputs"))
        vtk = root.find("Mesh/VTKMesh")
        vtk.set("logLevel", "2")
        vtk.set("partitionMethod", partition_method)
        vtk.set("partitionRefinement", str(partition_refinement))
        if mesh:
            vtk.set("file", str(mesh))
            # Collocation references are input main-point IDs and must stay exact.
            vtk.set("useGlobalIds", "1" if face_blocks else "0")
            vtk.set("mainBlockName", main_block)
            root.find("ElementRegions/CellElementRegion").set("cellBlocks", "{ * }")
            if face_blocks:
                vtk.set("faceBlocks", "{" + ",".join(face_blocks) + "}")
                for index, block in enumerate(face_blocks):
                    ET.SubElement(root.find("ElementRegions"), "SurfaceElementRegion",
                                  name=f"fracture{index}", faceBlock=block,
                                  materialList="{water,fractureRock}", defaultAperture="1e-4")
                root.find("Solvers/SinglePhaseFVM").set(
                    "targetRegions", "{region," + ",".join(f"fracture{i}" for i in range(len(face_blocks))) + "}")
                constitutive = root.find("Constitutive")
                ET.SubElement(constitutive, "CompressibleSolidParallelPlatesPermeability", name="fractureRock",
                              solidModelName="nullSolid", porosityModelName="fracturePorosity",
                              permeabilityModelName="fracturePermeability")
                ET.SubElement(constitutive, "PressurePorosity", name="fracturePorosity",
                              defaultReferencePorosity="1", referencePressure="0", compressibility="0")
                ET.SubElement(constitutive, "ParallelPlatesPermeability", name="fracturePermeability")
            fields = root.find("FieldSpecifications")
            for field in list(fields):
                if field.get("initialCondition") != "1":
                    fields.remove(field)
            ET.SubElement(fields, "SourceFlux", name="source", objectPath="ElementRegions/region",
                          scale="-1e-3", setNames="{ all }")
            root.remove(root.find("Geometry"))
        tree.write(path, encoding="utf-8", xml_declaration=True)
        paths[level] = path
    return paths


def summarize(records, directory):
    def median_phase(samples, name):
        values = [r.get("phase_seconds", {}).get(name) for r in samples]
        return statistics.median(values) if all(v is not None for v in values) else None

    def median_metric(samples, section, name):
        values = [r["refinement_metrics"][section][-1][name] for r in samples
                  if r.get("refinement_metrics", {}).get(section)]
        return statistics.median(values) if len(values) == len(samples) else None

    groups = {}
    for record in records:
        if record["returncode"] == 0:
            groups.setdefault((record["mode"], record["nodes"], record["ranks"], record["level"]), []).append(record)
    rows = []
    for key, samples in sorted(groups.items()):
        mode, nodes, ranks, level = key
        baseline_key = min((k for k in groups if k[0] == mode and k[3] == level), key=lambda k: k[2])
        baseline = statistics.median(r["maximum_rank_seconds"] for r in groups[baseline_key])
        seconds = statistics.median(r["maximum_rank_seconds"] for r in samples)
        refinement = median_phase(samples, "refinement")
        baseline_refinement = median_phase(groups[baseline_key], "refinement")
        refinement_efficiency = (baseline_refinement / refinement *
                                 (baseline_key[2] / ranks if mode == "strong" else 1)
                                 if refinement and baseline_refinement is not None else None)
        rows.append({"mode": mode, "nodes": nodes, "ranks": ranks, "level": level,
                     "repeats": len(samples), "median_application_seconds": seconds,
                     "min_application_seconds": min(r["maximum_rank_seconds"] for r in samples),
                     "max_application_seconds": max(r["maximum_rank_seconds"] for r in samples),
                     "baseline_ranks": baseline_key[2],
                     "application_efficiency": baseline / seconds * (baseline_key[2] / ranks if mode == "strong" else 1),
                     "median_refinement_seconds": refinement,
                     "refinement_efficiency": refinement_efficiency,
                     "median_refinement_reporting_seconds": median_phase(samples, "refinement_reporting"),
                     "median_fracture_association_seconds": median_phase(samples, "fracture_associations"),
                     "median_connectivity_seconds": median_phase(samples, "connectivity"),
                     "median_ghost_setup_seconds": median_phase(samples, "ghosts"),
                     "median_dof_setup_seconds": median_phase(samples, "dofs"),
                     "median_sparsity_setup_seconds": median_phase(samples, "sparsity"),
                     "median_system_setup_seconds": median_phase(samples, "system_setup"),
                     "median_linear_iterations": statistics.median(r["linear_iterations"] for r in samples),
                     "median_rank0_initialization_seconds": statistics.median(r["rank0_initialization_seconds"] for r in samples),
                     "peak_rank_rss_kib": max(r["peak_rank_rss_kib"] for r in samples),
                     "final_owned_cells": median_metric(samples, "levels", "ownedCells"),
                     "final_max_rank_owned_cells": median_metric(samples, "levels", "maxOwnedCells"),
                     "final_cell_imbalance": median_metric(samples, "levels", "cellImbalance"),
                     "final_unique_main_points": median_metric(samples, "levels", "uniqueMainPoints"),
                     "final_main_point_copies": median_metric(samples, "levels", "mainPointCopies"),
                     "final_shared_main_point_copies": median_metric(samples, "levels", "sharedMainPointCopies"),
                     "final_max_rank_payload_messages": median_metric(samples, "levels", "maxRankPayloadMessages"),
                     "final_max_rank_payload_bytes": median_metric(samples, "levels", "maxRankPayloadBytes"),
                     "final_max_rank_count_bytes": median_metric(samples, "levels", "maxRankCountBytes"),
                     "final_max_rank_neighbor_exchanges": median_metric(samples, "levels", "maxRankNeighborExchanges"),
                     "final_max_point_copies_bound": median_metric(samples, "forecasts", "maxPointCopiesBound"),
                     "final_max_field_bytes_bound": median_metric(samples, "forecasts", "maxFieldBytesBound"),
                     "final_max_vtk_bytes_bound": median_metric(samples, "forecasts", "maxVtkBytesBound"),
                     "final_max_refiner_peak_bytes_model": median_metric(samples, "forecasts", "maxRefinerPeakBytesModel"),
                     "final_max_geos_owned_connectivity_bytes_model": median_metric(samples, "forecasts", "maxGeosOwnedConnectivityBytesModel"),
                     "final_max_exchange_bytes_bound": median_metric(samples, "forecasts", "maxExchangeBytesBound"),
                     "final_max_geos_ghost_connectivity_bytes_model": median_metric(samples, "forecasts", "maxGeosGhostConnectivityBytesModel")})
    if rows:
        with (directory / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def refinement_metrics(text, level, ranks):
    """Parse and check production reports, retaining each fine level separately."""
    forecasts, levels = {}, {}
    forecast_keys = {"ownedCells", "maxOwnedCells", "maxPointCopiesBound", "maxFieldBytesBound",
                     "maxVtkBytesBound", "maxRefinerPeakBytesModel", "maxGeosOwnedConnectivityBytesModel",
                     "maxExchangeBytesBound", "maxGeosGhostConnectivityBytesModel"}
    level_keys = {"ownedCells", "maxOwnedCells", "meanOwnedCells", "uniqueMainPoints", "mainPointCopies",
                  "sharedMainPointCopies", "maxRankDirectoryExchanges", "maxRankNeighborExchanges",
                  "maxRankPayloadMessages", "maxRankPayloadBytes", "maxRankCountBytes"}
    floating = {"meanOwnedCells", "maxRefinerPeakBytesModel", "maxExchangeBytesBound",
                "maxGeosGhostConnectivityBytesModel"}
    for forecast, number, body in re.findall(r"Uniform refinement (forecast )?level (\d+): ([^\n]+)", text):
        target = forecasts if forecast else levels
        number = int(number)
        if number in target:
            raise RuntimeError(f"Duplicate refinement report for level {number}")
        values = {}
        for item in body.strip().split(", "):
            key, value = item.split("=", 1)
            values[key] = float(value) if key in floating else int(value)
            if not math.isfinite(values[key]) or values[key] < 0:
                raise RuntimeError(f"Invalid refinement metric {key}={value}")
        if values.keys() != (forecast_keys if forecast else level_keys):
            raise RuntimeError(f"Unexpected refinement report fields: {values.keys()}")
        target[number] = {"level": number, **values}
    if set(forecasts) != set(range(1, level + 1)) or set(levels) != (set(range(level + 1)) if level else set()):
        raise RuntimeError(f"Incomplete refinement reports for requested level {level}")
    for number, values in levels.items():
        mean = values["ownedCells"] / ranks
        if not mean or not math.isclose(values["meanOwnedCells"], mean, rel_tol=1e-12):
            raise RuntimeError(f"Invalid owned-cell mean at refinement level {number}")
        if not mean <= values["maxOwnedCells"] <= values["ownedCells"]:
            raise RuntimeError(f"Invalid owned-cell maximum at refinement level {number}")
        copies, shared, unique = (values[k] for k in ("mainPointCopies", "sharedMainPointCopies", "uniqueMainPoints"))
        if not 0 <= shared <= copies or not copies - shared <= unique <= copies:
            raise RuntimeError(f"Invalid shared-point counts at refinement level {number}")
        values["cellImbalance"] = values["maxOwnedCells"] / mean
        if number:
            forecast = forecasts[number]
            if any(forecast[k] != values[k] for k in ("ownedCells", "maxOwnedCells")):
                raise RuntimeError(f"Refined cells differ from forecast at level {number}")
            if values["maxRankDirectoryExchanges"]:
                raise RuntimeError(f"Fine refinement used a directory exchange at level {number}")
    return {"forecasts": [forecasts[k] for k in sorted(forecasts)],
            "levels": [levels[k] for k in sorted(levels)]}


def profile_rows(path, case):
    """Retain the hierarchy and the inclusive min/max/mean rank timings."""
    stack, rows = [], []
    for line in path.read_text().splitlines():
        match = re.fullmatch(r"( *)(.*?)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*", line)
        if not match:
            continue
        indent, name, minimum, maximum, mean, _ = match.groups()
        depth = len(indent) // 2
        stack = stack[:depth]
        stack.append(name.strip('"'))
        rows.append({"case": case, "depth": depth, "scope": stack[-1],
                     "scope_path": "/".join(stack), "seconds_min": float(minimum),
                     "seconds_max": float(maximum), "seconds_mean": float(mean)})
    if not rows:
        raise RuntimeError(f"No runtime-report timing rows: {path}")
    return rows


def profile_totals(rows, level):
    scopes = {"refinement": "uniformRefinement/total", "refinement_reporting": "uniformRefinement/statistics",
              "fracture_associations": "uniformRefinement/fractureAssociations",
              "connectivity": "geos::CellBlockManager::buildMaps",
              "ghosts": "geos::CommunicationTools::setupGhosts",
              "dofs": "geos/DOFSetup", "sparsity": "geos/sparsitySetup",
              "system_setup": "geos::PhysicsSolverBase::setupSystem"}
    totals = {}
    for name, scope in scopes.items():
        values = [r["seconds_max"] for r in rows if r["scope"] == scope]
        totals[name] = sum(values) if values else (0 if name == "fracture_associations" or
                                                  (name in ("refinement", "refinement_reporting") and level == 0) else None)
    return totals


def mesh_hashes(mesh):
    """Hash the container and referenced pieces, including nested VTM/PVTU files."""
    hashes, pending = {}, [mesh]
    while pending:
        path = pending.pop().resolve(strict=True)
        if str(path) in hashes:
            continue
        hasher = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                hasher.update(chunk)
        hashes[str(path)] = hasher.hexdigest()
        if path.suffix.lower() in (".vtm", ".pvtu"):
            for entry in ET.parse(path).iter():
                reference = entry.get("file") or entry.get("Source")
                if reference:
                    pending.append(path.parent / reference)
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path, help="new shared result directory")
    parser.add_argument("--launcher", default="mpiexec")
    parser.add_argument("--nodes", type=positive, nargs="+", default=[1, 2, 4])
    parser.add_argument("--ranks-per-node", type=positive, required=True)
    parser.add_argument("--strong-cells", type=positive, nargs=3, default=[32, 32, 32])
    parser.add_argument("--weak-cells", type=positive, nargs=3, default=[4, 4, 4])
    parser.add_argument("--levels", type=int, choices=[0, 1, 2], nargs="+", default=[0, 1, 2])
    parser.add_argument("--modes", choices=["strong", "weak"], nargs="+", default=["strong", "weak"])
    parser.add_argument("--mesh", type=Path, help="existing VTU or VTM, strong scaling only")
    parser.add_argument("--main-block", default="main", help="main volume block in a VTM")
    parser.add_argument("--face-blocks", nargs="+", default=[], help="VTM fracture blocks; requires input global IDs")
    parser.add_argument("--partition-method", choices=["parmetis", "ptscotch"], default="parmetis")
    parser.add_argument("--partition-refinement", type=int, default=0)
    parser.add_argument("--repeats", type=positive, default=3)
    parser.add_argument("--timeout", type=positive, default=1800)
    parser.add_argument("--timers", default="runtime-report,calc.inclusive,output=caliper-report.txt,max_column_width=200",
                        help="GEOS -t configuration; requires a Caliper-enabled build")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    output = args.output.resolve()
    launcher = shlex.split(args.launcher)
    slurm = bool(launcher) and Path(launcher[0]).name == "srun"
    if not launcher or not os.access(binary, os.X_OK):
        parser.error("provide an executable binary and a nonempty launcher")
    if args.partition_refinement < 0:
        parser.error("--partition-refinement must be nonnegative")
    if args.mesh and args.modes != ["strong"]:
        parser.error("--mesh requires --modes strong")
    if args.face_blocks and (not args.mesh or args.mesh.suffix.lower() != ".vtm"):
        parser.error("--face-blocks requires --mesh FILE.vtm")
    if len(set(args.face_blocks)) != len(args.face_blocks) or args.main_block in args.face_blocks:
        parser.error("face block names must be distinct and exclude the main block")
    if any(not re.fullmatch(r"[A-Za-z0-9_./*\[\]-]+", name) for name in [args.main_block, *args.face_blocks]):
        parser.error("block names must be GEOS group references without commas")
    if slurm and not args.dry_run and "SLURM_JOB_ID" not in os.environ:
        parser.error("use srun inside an existing allocation")
    if not Path("/usr/bin/time").is_file():
        parser.error("GNU /usr/bin/time is required on all compute nodes")
    mesh = args.mesh.resolve(strict=True) if args.mesh else None
    output.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parent.parent
    metadata = provenance(binary, repo, "boundary", "interfaces")
    for path in [Path(__file__).resolve(), repo / "benchmarks/verifyUniformRefinementInitialization.py",
                 *repo.glob("src/coreComponents/mesh/generators/VTKUniformRefinement*.[ch]pp"),
                 *repo.glob("src/coreComponents/mesh/generators/VTKMeshGenerator*.[ch]pp"),
                 *repo.glob("src/coreComponents/mesh/generators/VTKUtilities*.[ch]pp"),
                 *repo.glob("src/coreComponents/mesh/WellElementSubRegion*.[ch]pp"),
                 *repo.glob("src/coreComponents/mesh/FaceElementSubRegion*.[ch]pp"),
                 *repo.glob("src/coreComponents/mesh/mpiCommunications/SpatialPartition*.[ch]pp"),
                 repo / "src/coreComponents/mesh/MeshLevel.cpp",
                 repo / "src/coreComponents/mesh/DomainPartition.cpp",
                 repo / "src/coreComponents/mesh/ElementRegionManager.cpp",
                 repo / "src/coreComponents/physicsSolvers/PhysicsSolverBase.cpp"]:
        metadata["source_sha256"][str(path.relative_to(repo))] = hashlib.sha256(path.read_bytes()).hexdigest()
    for name in ("libvtkUniformRefinement.so", "libphysicsSolversBase.so"):
        path = binary.parent.parent / "lib" / name
        if path.is_file():
            metadata["component_library_sha256"][str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    if mesh:
        metadata["input_mesh_files_sha256"] = mesh_hashes(mesh)
        metadata["input_mesh_sha256"] = metadata["input_mesh_files_sha256"][str(mesh)]
    metadata.update(component_only=False, measured_path="full initialization and one FVM step, no checkpoint I/O",
                    arguments=vars(args) | {"binary": str(binary), "output": str(output), "mesh": str(mesh) if mesh else None})
    (output / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    records = []
    phases = []
    for mode in dict.fromkeys(args.modes):
        for nodes in sorted(set(args.nodes)):
            ranks = nodes * args.ranks_per_node
            grid = process_grid(ranks)
            cells = tuple(a * b for a, b in zip(args.weak_cells, grid)) if mode == "weak" else tuple(args.strong_cells)
            domain = tuple(a * b for a, b in zip((4, 2, 2), grid)) if mode == "weak" else (4, 2, 2)
            inputs = make_case_inputs(output / f"inputs-{mode}-n{nodes}", cells, domain, mesh,
                                      args.partition_method, args.partition_refinement,
                                      args.main_block, args.face_blocks)
            for level in sorted(set(args.levels)):
                for repeat in range(args.repeats):
                    label = f"{mode}-n{nodes}-p{ranks}-l{level}-r{repeat}"
                    directory = output / label
                    directory.mkdir()
                    geos = [str(binary), "-i", str(inputs[level]), "-o", str(directory)]
                    if args.timers:
                        geos += ["-t", args.timers]
                    placement = ([f"--nodes={nodes}", f"--ntasks={ranks}", f"--ntasks-per-node={args.ranks_per_node}"]
                                 if slurm else ["-n", str(ranks)])
                    command = launcher + placement + [sys.executable, str(Path(__file__).resolve()),
                                                       "--rank-wrapper", str(directory), *geos]
                    record = {"label": label, "mode": mode, "nodes": nodes, "ranks": ranks,
                              "level": level, "repeat": repeat, "global_coarse_cells": None if mesh else math.prod(cells),
                              "physical_domain": None if mesh else domain,
                              "command": command, "returncode": "dry-run"}
                    print(shlex.join(command), flush=True)
                    if not args.dry_run:
                        with (directory / "geos.log").open("w") as log:
                            process = subprocess.Popen(command, cwd=directory, stdout=log, stderr=log, start_new_session=True)
                            try:
                                record["returncode"] = process.wait(timeout=args.timeout)
                            except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                                record["returncode"] = "timeout" if isinstance(error, subprocess.TimeoutExpired) else "interrupted"
                                try:
                                    os.killpg(process.pid, signal.SIGTERM)
                                except ProcessLookupError:
                                    pass
                                try:
                                    process.wait(timeout=10)
                                except subprocess.TimeoutExpired:
                                    os.killpg(process.pid, signal.SIGKILL)
                                    process.wait()
                        if record["returncode"] == 0:
                            timings = [path.read_text().strip().split(",") for path in sorted(directory.glob("rank-*.time"))]
                            if len(timings) != ranks or any(len(t) != 3 or int(t[2]) != 0 for t in timings):
                                raise RuntimeError(f"Incomplete per-rank measurements: {directory}")
                            record["maximum_rank_seconds"] = max(float(t[0]) for t in timings)
                            record["peak_rank_rss_kib"] = max(int(t[1]) for t in timings)
                            placements = [json.loads(path.read_text()) for path in sorted(directory.glob("rank-*.json"))]
                            hosts = Counter(p["host"] for p in placements)
                            record["actual_ranks_per_host"] = dict(hosts)
                            if slurm and (len(hosts) != nodes or set(hosts.values()) != {args.ranks_per_node}):
                                raise RuntimeError(f"Unexpected Slurm placement: {hosts}")
                            text = (directory / "geos.log").read_text()
                            record["refinement_metrics"] = refinement_metrics(text, level, ranks)
                            record["linear_iterations"] = sum(int(value) for value in re.findall(
                                r"Linear solve:\s+\( iter, res \) = \(\s*(\d+),", text))
                            if record["linear_iterations"] == 0:
                                raise RuntimeError(f"Expected a nonzero FVM solve: {directory / 'geos.log'}")
                            record["geos_time_lines"] = re.findall(r"^(?:total|initialization|run) time.*$", text, re.MULTILINE)
                            record["mesh_count_lines"] = re.findall(r"^.*total elements.*$", text, re.MULTILINE)
                            initialization = re.search(r"^initialization time\s+.*\(([0-9.eE+-]+) s\)", text, re.MULTILINE)
                            if not initialization:
                                raise RuntimeError(f"Missing GEOS initialization time: {directory}")
                            record["rank0_initialization_seconds"] = float(initialization.group(1))
                            if "output=caliper-report.txt" in args.timers:
                                report = directory / "caliper-report.txt"
                                if not report.is_file() or "geos::ProblemManager::problemSetup" not in report.read_text():
                                    raise RuntimeError(f"Missing Caliper initialization profile: {directory}")
                                if level and "uniformRefinement/level" not in report.read_text():
                                    raise RuntimeError(f"Missing Caliper refinement profile: {directory}")
                                current_phases = profile_rows(report, label)
                                record["phase_seconds"] = profile_totals(current_phases, level)
                                if level and args.face_blocks and not record["phase_seconds"]["fracture_associations"]:
                                    raise RuntimeError(f"Missing fracture-association timing: {directory}")
                                phases.extend(current_phases)
                                with (output / "phases.csv").open("w", newline="") as stream:
                                    writer = csv.DictWriter(stream, fieldnames=list(phases[0]))
                                    writer.writeheader()
                                    writer.writerows(phases)
                            if not record["geos_time_lines"] or not record["mesh_count_lines"]:
                                raise RuntimeError(f"Missing completed GEOS initialization/count report: {directory}")
                    records.append(record)
                    (output / "runs.json").write_text(json.dumps(records, indent=2) + "\n")
                    summarize(records, output)
                    if record["returncode"] not in (0, "dry-run"):
                        raise SystemExit(f"GEOS failed ({record['returncode']}): {directory / 'geos.log'}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--rank-wrapper":
        rank_wrapper(sys.argv[2:])
    else:
        main()
