#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

"""Build a vendor-independent, advisory pull-request coverage report."""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import pathlib
import re
import subprocess
import sys
import tempfile
from typing import Iterable


SCRIPT_DIRECTORY = pathlib.Path( __file__ ).resolve().parent
if str( SCRIPT_DIRECTORY ) not in sys.path:
    sys.path.insert( 0, str( SCRIPT_DIRECTORY ) )
import check_coverage_thresholds as coverage_contract


REPORT_SCHEMA_VERSION = 1
SUPPORTED_SUMMARY_SCHEMA_VERSIONS = ( 2, 3 )
PRODUCTION_SCOPE = "src/coreComponents"
CANONICAL_METRICS = ( "regions", "functions", "lines", "branches" )
DISPLAY_METRICS = ( "lines", "functions", "branches", "regions" )
SHA_PATTERN = re.compile( r"^[0-9a-fA-F]{40}$" )
HUNK_PATTERN = re.compile(
    rb"^@@ -[0-9]+(?:,[0-9]+)? \+([0-9]+)(?:,([0-9]+))? @@"
)
MAX_IMPACTED_FILES = 10
MAX_IMPACTED_LOCATIONS = 10
MAX_FUNCTION_NAMES = 3


def require_integer( value, description: str, minimum: int = 0 ) -> int:
    if isinstance( value, bool ) or not isinstance( value, int ) or value < minimum:
        raise ValueError( f"{description} must be an integer >= {minimum}" )
    return value


def validate_sha( value: str, description: str ) -> str:
    if not isinstance( value, str ) or SHA_PATTERN.fullmatch( value ) is None:
        raise ValueError( f"{description} must be a full 40-character Git SHA" )
    return value.lower()


def read_json( path: pathlib.Path, description: str ) -> dict:
    try:
        with path.open( encoding="utf-8" ) as stream:
            document = json.load( stream )
    except ( OSError, UnicodeError, json.JSONDecodeError ) as error:
        raise ValueError( f"cannot read {description} {path}: {error}" ) from error
    if not isinstance( document, dict ):
        raise ValueError( f"{description} must contain a JSON object" )
    return document


def validate_metric( raw_metric: object, description: str ) -> dict:
    if not isinstance( raw_metric, dict ):
        raise ValueError( f"{description} must be an object" )
    covered = require_integer( raw_metric.get( "covered" ), f"{description}.covered" )
    total = require_integer( raw_metric.get( "total" ), f"{description}.total", 1 )
    not_covered = require_integer(
        raw_metric.get( "not_covered" ), f"{description}.not_covered"
    )
    if covered > total or not_covered != total - covered:
        raise ValueError( f"{description} counts are inconsistent" )
    percent = raw_metric.get( "percent" )
    if isinstance( percent, bool ) or not isinstance( percent, ( int, float ) ):
        raise ValueError( f"{description}.percent must be numeric" )
    expected_percent = round( 100.0 * covered / total, 6 )
    if not math.isfinite( float( percent ) ) or abs(
        float( percent ) - expected_percent
    ) > 1.0e-6:
        raise ValueError( f"{description}.percent is inconsistent with its counts" )
    return {
        "covered": covered,
        "not_covered": not_covered,
        "total": total,
        "percent": expected_percent,
    }


def validate_summary( document: dict, description: str ) -> dict:
    schema_version = document.get( "schema_version" )
    if (
        isinstance( schema_version, bool )
        or not isinstance( schema_version, int )
        or schema_version not in SUPPORTED_SUMMARY_SCHEMA_VERSIONS
    ):
        supported = ", ".join( str( value ) for value in SUPPORTED_SUMMARY_SCHEMA_VERSIONS )
        raise ValueError(
            f"{description}.schema_version must be one of the integers {supported}"
        )
    if schema_version == coverage_contract.COVERAGE_SUMMARY_SCHEMA_VERSION:
        try:
            validated = coverage_contract.validate_summary( document )
        except ValueError as error:
            raise ValueError( f"{description} is invalid: {error}" ) from error
        return {
            "schema_version": schema_version,
            "scope": validated["scope"],
            "excluded_regex": validated["excluded_regex"],
            "tool": validated["tool"],
            "metrics": validated["metrics"],
            "measurement": validated["measurement"],
            "per_file_metrics": validated["per_file_metrics"],
        }
    scope = document.get( "scope" )
    excluded_regex = document.get( "excluded_regex" )
    if not isinstance( scope, str ) or not scope:
        raise ValueError( f"{description}.scope must be a nonempty string" )
    if not isinstance( excluded_regex, str ) or not excluded_regex:
        raise ValueError( f"{description}.excluded_regex must be a nonempty string" )
    tool = document.get( "tool" )
    if not isinstance( tool, dict ) or tool.get( "name" ) != "llvm-cov":
        raise ValueError( f"{description}.tool.name must be llvm-cov" )
    tool_major = require_integer( tool.get( "major" ), f"{description}.tool.major", 1 )
    raw_metrics = document.get( "metrics" )
    if not isinstance( raw_metrics, dict ):
        raise ValueError( f"{description}.metrics must be an object" )
    metrics = {
        name: validate_metric( raw_metrics.get( name ), f"{description}.metrics.{name}" )
        for name in CANONICAL_METRICS
    }

    raw_measurement = document.get( "measurement" )
    measurement = None
    if raw_measurement is not None:
        raise ValueError( f"{description}.measurement is not valid in schema 2" )

    return {
        "schema_version": schema_version,
        "scope": scope,
        "excluded_regex": excluded_regex,
        "tool": { "name": "llvm-cov", "major": tool_major },
        "metrics": metrics,
        "measurement": measurement,
        "per_file_metrics": None,
    }


def checked_input_file( path: pathlib.Path, description: str ) -> pathlib.Path:
    if path.is_symlink():
        raise ValueError( f"{description} must not be a symbolic link: {path}" )
    try:
        resolved = path.resolve( strict=True )
    except OSError as error:
        raise ValueError( f"cannot resolve {description} {path}: {error}" ) from error
    if not resolved.is_file():
        raise ValueError( f"{description} is not a regular file: {path}" )
    return resolved


def run_git( repository: pathlib.Path, arguments: Iterable[str] ) -> bytes:
    command = [
        "git",
        "-c",
        f"safe.directory={repository}",
        "-C",
        str( repository ),
        *arguments,
    ]
    try:
        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as error:
        raise ValueError( f"cannot execute Git: {error}" ) from error
    if result.returncode != 0:
        message = result.stderr.decode( "utf-8", errors="replace" ).strip()
        raise ValueError( f"Git command failed: {message or 'unknown error'}" )
    return result.stdout


def validate_repository( repository: pathlib.Path ) -> pathlib.Path:
    try:
        resolved = repository.resolve( strict=True )
    except OSError as error:
        raise ValueError( f"cannot resolve repository {repository}: {error}" ) from error
    if not resolved.is_dir():
        raise ValueError( f"repository is not a directory: {resolved}" )
    top_level = run_git( resolved, ( "rev-parse", "--show-toplevel" ) )
    try:
        git_root = pathlib.Path( top_level.decode( "utf-8" ).strip() ).resolve(
            strict=True
        )
    except ( UnicodeError, OSError ) as error:
        raise ValueError( f"cannot resolve Git repository root: {error}" ) from error
    if git_root != resolved:
        raise ValueError( f"--repository must name the Git top level: {git_root}" )
    return resolved


def validate_commit( repository: pathlib.Path, sha: str, description: str ) -> str:
    normalized = validate_sha( sha, description )
    resolved = run_git(
        repository, ( "rev-parse", "--verify", f"{normalized}^{{commit}}" )
    ).decode( "ascii", errors="strict" ).strip()
    if resolved.lower() != normalized:
        raise ValueError( f"{description} does not resolve to the requested commit" )
    return normalized


def commit_tree( repository: pathlib.Path, sha: str ) -> str:
    tree = run_git(
        repository, ( "rev-parse", "--verify", f"{sha}^{{tree}}" )
    ).decode( "ascii", errors="strict" ).strip()
    return validate_sha( tree, f"tree for {sha}" )


def decode_git_path( value: bytes ) -> str:
    try:
        path = value.decode( "utf-8", errors="strict" )
    except UnicodeError as error:
        raise ValueError( "Git paths must be valid UTF-8" ) from error
    normalized = pathlib.PurePosixPath( path )
    if (
        not path
        or normalized.is_absolute()
        or ".." in normalized.parts
        or normalized.as_posix() != path
    ):
        raise ValueError( f"Git returned an unsafe or non-normalized path: {path!r}" )
    return path


def changed_line_ranges(
    repository: pathlib.Path, base_sha: str, head_sha: str
) -> list[dict]:
    raw_names = run_git(
        repository,
        (
            "diff",
            "--name-status",
            "-z",
            "--find-renames",
            "--find-copies-harder",
            "--diff-filter=ACMR",
            base_sha,
            head_sha,
            "--",
            PRODUCTION_SCOPE,
        ),
    )
    fields = raw_names.split( b"\0" )
    if fields and fields[-1] == b"":
        fields.pop()

    records = []
    index = 0
    scope_prefix = f"{PRODUCTION_SCOPE}/"
    while index < len( fields ):
        try:
            status_text = fields[index].decode( "ascii", errors="strict" )
        except UnicodeError as error:
            raise ValueError( "Git returned a non-ASCII change status" ) from error
        index += 1
        if not status_text or status_text[0] not in "ACMR":
            raise ValueError( f"Git returned an unexpected change status: {status_text!r}" )
        status = status_text[0]
        if status in "RC":
            if index + 1 >= len( fields ):
                raise ValueError( "Git returned a truncated rename/copy record" )
            old_path = decode_git_path( fields[index] )
            path = decode_git_path( fields[index + 1] )
            index += 2
            diff_paths = ( old_path, path )
        else:
            if index >= len( fields ):
                raise ValueError( "Git returned a truncated change record" )
            path = decode_git_path( fields[index] )
            index += 1
            old_path = None
            diff_paths = ( path, )

        # A rename out of the production tree can be selected by its old path.
        # It has no head lines in scope and therefore is not a patch denominator.
        if not path.startswith( scope_prefix ):
            continue

        patch = run_git(
            repository,
            (
                "diff",
                "--unified=0",
                "--no-color",
                "--no-ext-diff",
                "--find-renames",
                "--find-copies-harder",
                base_sha,
                head_sha,
                "--",
                *diff_paths,
            ),
        )
        lines = set()
        for patch_line in patch.splitlines():
            match = HUNK_PATTERN.match( patch_line )
            if match is None:
                continue
            start = int( match.group( 1 ) )
            count = int( match.group( 2 ) or b"1" )
            if count:
                lines.update( range( start, start + count ) )
        records.append(
            {
                "status": status,
                "path": path,
                "old_path": old_path,
                "lines": frozenset( lines ),
            }
        )

    records.sort( key=lambda record: record["path"] )
    seen_paths = set()
    for record in records:
        if record["path"] in seen_paths:
            raise ValueError( f"duplicate changed path from Git: {record['path']}" )
        seen_paths.add( record["path"] )
    return records


def normalize_lcov_source( source: str, repository: pathlib.Path ) -> str | None:
    candidate = pathlib.Path( source )
    if not candidate.is_absolute():
        candidate = repository / candidate
    try:
        resolved = candidate.resolve( strict=False )
        relative = resolved.relative_to( repository )
    except ( OSError, ValueError ):
        return None
    path = relative.as_posix()
    if not path.startswith( f"{PRODUCTION_SCOPE}/" ):
        return None
    return path


def parse_nonnegative_integer( value: str, description: str, minimum: int = 0 ) -> int:
    try:
        parsed = int( value, 10 )
    except ValueError as error:
        raise ValueError( f"{description} must be an integer" ) from error
    if parsed < minimum:
        raise ValueError( f"{description} must be >= {minimum}" )
    return parsed


def parse_native_lcov(
    native_info: pathlib.Path,
    repository: pathlib.Path,
    changed_files: list[dict],
) -> dict:
    changed_lines = {
        record["path"]: set( record["lines"] ) for record in changed_files
    }
    line_counts: dict[tuple[str, int], int] = {}
    branch_outcomes: dict[tuple[str, int, int, int], bool] = {}
    function_names: dict[tuple[str, int], set[str]] = {}
    relevant_function_names: dict[tuple[str, str], set[int]] = {}
    function_counts: dict[tuple[str, str], int] = {}
    current_path = None

    try:
        stream = native_info.open( encoding="utf-8", errors="strict" )
    except OSError as error:
        raise ValueError( f"cannot read native LCOV input {native_info}: {error}" ) from error
    try:
        with stream:
            for line_number, raw_line in enumerate( stream, 1 ):
                line = raw_line.rstrip( "\r\n" )
                if line.startswith( "SF:" ):
                    current_path = normalize_lcov_source( line[3:], repository )
                    if current_path not in changed_lines:
                        current_path = None
                    continue
                if line == "end_of_record":
                    current_path = None
                    continue
                if current_path is None:
                    continue

                if line.startswith( "DA:" ):
                    fields = line[3:].split( ",", 2 )
                    if len( fields ) < 2:
                        raise ValueError( f"malformed DA record at LCOV line {line_number}" )
                    source_line = parse_nonnegative_integer(
                        fields[0], f"DA source line at LCOV line {line_number}", 1
                    )
                    count = parse_nonnegative_integer(
                        fields[1], f"DA count at LCOV line {line_number}"
                    )
                    if source_line in changed_lines[current_path]:
                        key = ( current_path, source_line )
                        line_counts[key] = max( line_counts.get( key, 0 ), count )
                    continue

                if line.startswith( "FN:" ):
                    fields = line[3:].split( ",", 1 )
                    if len( fields ) != 2 or not fields[1]:
                        raise ValueError( f"malformed FN record at LCOV line {line_number}" )
                    source_line = parse_nonnegative_integer(
                        fields[0], f"FN source line at LCOV line {line_number}", 1
                    )
                    if source_line in changed_lines[current_path]:
                        name = fields[1]
                        function_names.setdefault(
                            ( current_path, source_line ), set()
                        ).add( name )
                        relevant_function_names.setdefault(
                            ( current_path, name ), set()
                        ).add( source_line )
                    continue

                if line.startswith( "FNDA:" ):
                    fields = line[5:].split( ",", 1 )
                    if len( fields ) != 2 or not fields[1]:
                        raise ValueError( f"malformed FNDA record at LCOV line {line_number}" )
                    name = fields[1]
                    if ( current_path, name ) not in relevant_function_names:
                        continue
                    count = parse_nonnegative_integer(
                        fields[0], f"FNDA count at LCOV line {line_number}"
                    )
                    key = ( current_path, name )
                    function_counts[key] = max( function_counts.get( key, 0 ), count )
                    continue

                if line.startswith( "BRDA:" ):
                    fields = line[5:].split( ",", 3 )
                    if len( fields ) != 4:
                        raise ValueError( f"malformed BRDA record at LCOV line {line_number}" )
                    source_line = parse_nonnegative_integer(
                        fields[0], f"BRDA source line at LCOV line {line_number}", 1
                    )
                    if source_line not in changed_lines[current_path]:
                        continue
                    pair_id = parse_nonnegative_integer(
                        fields[1], f"BRDA pair id at LCOV line {line_number}"
                    )
                    branch_id = parse_nonnegative_integer(
                        fields[2], f"BRDA branch id at LCOV line {line_number}"
                    )
                    if fields[3] == "-":
                        covered = False
                    else:
                        covered = parse_nonnegative_integer(
                            fields[3], f"BRDA count at LCOV line {line_number}"
                        ) > 0
                    key = ( current_path, source_line, pair_id, branch_id )
                    branch_outcomes[key] = branch_outcomes.get( key, False ) or covered
    except UnicodeError as error:
        raise ValueError( f"native LCOV input must be valid UTF-8: {error}" ) from error

    function_sites = {}
    for site, names in function_names.items():
        covered = any(
            function_counts.get( ( site[0], name ), 0 ) > 0 for name in names
        )
        function_sites[site] = {
            "covered": covered,
            "names": sorted( names ),
        }
    return {
        "line_counts": line_counts,
        "branch_outcomes": branch_outcomes,
        "function_sites": function_sites,
    }


def patch_metric( covered: int, total: int ) -> dict:
    return {
        "covered": covered,
        "not_covered": total - covered,
        "total": total,
        "percent": round( 100.0 * covered / total, 6 ) if total else None,
    }


def location_hint(
    line_missed: bool,
    missed_branches: int,
    covered_branches: int,
    function_missed: bool,
) -> str:
    if line_missed and function_missed:
        if missed_branches:
            return "Call this function in a coverage smoke, then exercise its decision outcomes."
        return "Add a coverage smoke that calls this function entry."
    if line_missed:
        if missed_branches:
            return "Reach this executable line, then exercise its decision outcomes."
        return "Add a coverage smoke that reaches this executable line."
    if missed_branches:
        if covered_branches:
            return "Add cases that exercise the remaining decision outcome(s)."
        return "Add cases that take each unobserved decision outcome."
    if function_missed:
        return "Add a coverage smoke that calls this function entry."
    return "Add a focused coverage smoke for this location."


def build_patch_report( changed_files: list[dict], lcov: dict ) -> dict:
    line_counts = lcov["line_counts"]
    branch_outcomes = lcov["branch_outcomes"]
    function_sites = lcov["function_sites"]
    file_results = []
    locations = []

    for changed_file in changed_files:
        path = changed_file["path"]
        changed = set( changed_file["lines"] )
        file_lines = {
            line: count
            for ( source_path, line ), count in line_counts.items()
            if source_path == path
        }
        file_branches = {
            key: covered
            for key, covered in branch_outcomes.items()
            if key[0] == path
        }
        file_functions = {
            line: details
            for ( source_path, line ), details in function_sites.items()
            if source_path == path
        }
        line_metric = patch_metric(
            sum( count > 0 for count in file_lines.values() ), len( file_lines )
        )
        function_metric = patch_metric(
            sum( details["covered"] for details in file_functions.values() ),
            len( file_functions ),
        )
        branch_metric = patch_metric(
            sum( file_branches.values() ), len( file_branches )
        )
        file_results.append(
            {
                "path": path,
                "old_path": changed_file["old_path"],
                "status": changed_file["status"],
                "changed_lines": len( changed ),
                "non_instrumented_changed_lines": len( changed - set( file_lines ) ),
                "metrics": {
                    "executable_lines": line_metric,
                    "function_start_sites": function_metric,
                    "native_branch_outcomes": branch_metric,
                },
            }
        )

        candidate_lines = set( file_lines ) | set( file_functions )
        candidate_lines.update( key[1] for key in file_branches )
        for line in candidate_lines:
            executable = line in file_lines
            line_covered = file_lines.get( line, 0 ) > 0 if executable else None
            branch_values = [
                covered for key, covered in file_branches.items() if key[1] == line
            ]
            covered_branches = sum( branch_values )
            missed_branches = len( branch_values ) - covered_branches
            function = file_functions.get( line )
            function_missed = function is not None and not function["covered"]
            line_missed = executable and not line_covered
            if not ( line_missed or missed_branches or function_missed ):
                continue
            names = function["names"] if function is not None else []
            locations.append(
                {
                    "path": path,
                    "line": line,
                    "executable_line": executable,
                    "line_covered": line_covered,
                    "native_branch_outcomes": patch_metric(
                        covered_branches, len( branch_values )
                    ),
                    "function_start_covered": (
                        function["covered"] if function is not None else None
                    ),
                    "function_names": names[:MAX_FUNCTION_NAMES],
                    "additional_function_names": max(
                        0, len( names ) - MAX_FUNCTION_NAMES
                    ),
                    "hint": location_hint(
                        line_missed,
                        missed_branches,
                        covered_branches,
                        function_missed,
                    ),
                }
            )

    file_results.sort(
        key=lambda result: (
            -result["metrics"]["native_branch_outcomes"]["not_covered"],
            -result["metrics"]["executable_lines"]["not_covered"],
            -result["metrics"]["function_start_sites"]["not_covered"],
            result["path"],
        )
    )
    locations.sort(
        key=lambda location: (
            -location["native_branch_outcomes"]["not_covered"],
            -( location["executable_line"] and not location["line_covered"] ),
            -( location["function_start_covered"] is False ),
            location["path"],
            location["line"],
        )
    )

    executable_total = len( line_counts )
    function_total = len( function_sites )
    branch_total = len( branch_outcomes )
    all_changed_lines = sum( len( record["lines"] ) for record in changed_files )
    return {
        "changed_files": len( changed_files ),
        "changed_lines": all_changed_lines,
        "non_instrumented_changed_lines": all_changed_lines - executable_total,
        "metrics": {
            "executable_lines": patch_metric(
                sum( count > 0 for count in line_counts.values() ), executable_total
            ),
            "function_start_sites": patch_metric(
                sum( details["covered"] for details in function_sites.values() ),
                function_total,
            ),
            "native_branch_outcomes": patch_metric(
                sum( branch_outcomes.values() ), branch_total
            ),
        },
        "top_impacted_files": file_results[:MAX_IMPACTED_FILES],
        "top_impacted_locations": locations[:MAX_IMPACTED_LOCATIONS],
    }


def unavailable_comparison( reason: str, candidate: dict ) -> dict:
    return {
        "comparable": False,
        "reason": reason,
        "baseline": None,
        "candidate": candidate["metrics"],
        "metrics": {},
        "changed_file_metrics": [],
        "provenance": None,
    }


def build_aggregate_comparison(
    candidate: dict,
    baseline: dict | None,
    baseline_error: str | None,
    base_sha: str,
    head_sha: str,
    base_tree: str,
    head_tree: str,
    changed_files: list[dict],
) -> dict:
    if baseline_error is not None:
        return unavailable_comparison( baseline_error, candidate )
    if baseline is None:
        return unavailable_comparison(
            "No exact-base coverage summary was provided.", candidate
        )
    if baseline["schema_version"] != 3 or candidate["schema_version"] != 3:
        return unavailable_comparison(
            "N/A: exact comparison requires schema-3 measurement provenance "
            "on both summaries.",
            candidate,
        )

    compatibility_fields = (
        ( "scope", "production scopes differ" ),
        ( "excluded_regex", "source exclusion contracts differ" ),
    )
    for field, reason in compatibility_fields:
        if baseline[field] != candidate[field]:
            return unavailable_comparison( f"N/A: {reason}.", candidate )
    if baseline["tool"] != candidate["tool"]:
        return unavailable_comparison( "N/A: LLVM coverage tools differ.", candidate )

    baseline_measurement = baseline["measurement"]
    candidate_measurement = candidate["measurement"]
    if baseline_measurement is None or candidate_measurement is None:
        raise ValueError( "schema-3 coverage summaries are missing validated provenance" )
    if baseline_measurement["commit_sha"] != base_sha:
        return unavailable_comparison(
            "N/A: baseline measurement commit does not match --base-sha.", candidate
        )
    if candidate_measurement["commit_sha"] != head_sha:
        return unavailable_comparison(
            "N/A: candidate measurement commit does not match --head-sha.", candidate
        )
    if baseline_measurement["tree_sha"] != base_tree:
        return unavailable_comparison(
            "N/A: baseline measurement tree does not match the base commit.", candidate
        )
    if candidate_measurement["tree_sha"] != head_tree:
        return unavailable_comparison(
            "N/A: candidate measurement tree does not match the head commit.", candidate
        )
    if baseline_measurement["contract_id"] != candidate_measurement["contract_id"]:
        return unavailable_comparison(
            "N/A: coverage measurement contract IDs differ.", candidate
        )
    if (
        baseline_measurement["contract_fingerprint"]
        != candidate_measurement["contract_fingerprint"]
    ):
        return unavailable_comparison(
            "N/A: coverage measurement contract fingerprints differ.", candidate
        )

    comparisons = {}
    for name in CANONICAL_METRICS:
        base_metric = baseline["metrics"][name]
        head_metric = candidate["metrics"][name]
        comparisons[name] = {
            "baseline": base_metric,
            "candidate": head_metric,
            "covered_delta": head_metric["covered"] - base_metric["covered"],
            "not_covered_delta": (
                head_metric["not_covered"] - base_metric["not_covered"]
            ),
            "total_delta": head_metric["total"] - base_metric["total"],
            "percentage_point_delta": round(
                head_metric["percent"] - base_metric["percent"], 6
            ),
        }
    baseline_files = {
        source["path"]: source["metrics"]
        for source in baseline["per_file_metrics"]
    }
    candidate_files = {
        source["path"]: source["metrics"]
        for source in candidate["per_file_metrics"]
    }
    empty_metric = { "covered": 0, "not_covered": 0, "total": 0, "percent": 0.0 }
    changed_file_metrics = []
    for changed_file in changed_files:
        path = changed_file["path"]
        baseline_path = changed_file["old_path"] or path
        baseline_metrics = baseline_files.get( baseline_path, {} )
        candidate_metrics = candidate_files.get( path, {} )
        metric_changes = {}
        for name in CANONICAL_METRICS:
            base_metric = baseline_metrics.get( name, empty_metric )
            head_metric = candidate_metrics.get( name, empty_metric )
            metric_changes[name] = {
                "baseline": base_metric,
                "candidate": head_metric,
                "covered_delta": head_metric["covered"] - base_metric["covered"],
                "not_covered_delta": (
                    head_metric["not_covered"] - base_metric["not_covered"]
                ),
                "total_delta": head_metric["total"] - base_metric["total"],
            }
        changed_file_metrics.append(
            {
                "path": path,
                "baseline_path": baseline_path,
                "metrics": metric_changes,
            }
        )
    changed_file_metrics.sort(
        key=lambda source: (
            -source["metrics"]["branches"]["not_covered_delta"],
            -source["metrics"]["lines"]["not_covered_delta"],
            source["path"],
        )
    )
    return {
        "comparable": True,
        "reason": "Exact commit, tree, and coverage measurement contract match.",
        "baseline": baseline["metrics"],
        "candidate": candidate["metrics"],
        "metrics": comparisons,
        "changed_file_metrics": changed_file_metrics[:MAX_IMPACTED_FILES],
        "provenance": {
            "baseline_commit_sha": base_sha,
            "candidate_commit_sha": head_sha,
            "contract_id": candidate_measurement["contract_id"],
            "contract_fingerprint": candidate_measurement["contract_fingerprint"],
        },
    }


def markdown_escape( value: str ) -> str:
    escaped = html.escape( value, quote=True )
    return (
        escaped.replace( "|", "&#124;" )
        .replace( "\r", "&#13;" )
        .replace( "\n", "&#10;" )
    )


def percent_text( metric: dict ) -> str:
    if metric["total"] == 0:
        return "N/A"
    return f"{100.0 * metric['covered'] / metric['total']:.2f}%"


def compact_metric_text( metric: dict ) -> str:
    if metric["total"] == 0:
        return "N/A"
    return (
        f"{metric['covered']:,}/{metric['total']:,} "
        f"({percent_text( metric )})"
    )


def signed_integer( value: int ) -> str:
    return f"{value:+,}"


def build_markdown( report: dict ) -> str:
    aggregate = report["aggregate_comparison"]
    patch = report["patch"]
    display_names = {
        "regions": "Regions",
        "functions": "Functions",
        "lines": "Lines",
        "branches": "Canonical branches",
    }
    rows = [
        "### Pull-request coverage",
        "",
        f"Effective source diff: `{report['base_sha'][:12]}` → "
        f"`{report['head_sha'][:12]}` · Production scope: "
        f"`{markdown_escape( report['scope'] )}`",
        "",
        "#### Project coverage change",
        "",
    ]
    if aggregate["comparable"]:
        rows.extend(
            (
                "| Metric | Exact base | Candidate | Coverage change | "
                "Covered Δ | Missed Δ | Total Δ |",
                "|---|---:|---:|---:|---:|---:|---:|",
            )
        )
        for name in DISPLAY_METRICS:
            comparison = aggregate["metrics"][name]
            rows.append(
                f"| {display_names[name]} | "
                f"{compact_metric_text( comparison['baseline'] )} | "
                f"{compact_metric_text( comparison['candidate'] )} | "
                f"{comparison['percentage_point_delta']:+.2f} pp | "
                f"{signed_integer( comparison['covered_delta'] )} | "
                f"{signed_integer( comparison['not_covered_delta'] )} | "
                f"{signed_integer( comparison['total_delta'] )} |"
            )
        rows.extend(
            (
                "",
                "Comparison is exact because commit, tree, and measurement-contract "
                "provenance match.",
            )
        )
        changed_file_metrics = aggregate["changed_file_metrics"]
        if changed_file_metrics:
            rows.extend(
                (
                    "",
                    "##### Canonical movement in changed files",
                    "",
                    "| Source file | Base branches | Candidate branches | "
                    "Missed branch Δ | Missed line Δ |",
                    "|---|---:|---:|---:|---:|",
                )
            )
            for source in changed_file_metrics:
                branch_change = source["metrics"]["branches"]
                line_change = source["metrics"]["lines"]
                display_path = source["path"]
                if source["baseline_path"] != source["path"]:
                    display_path = f"{source['baseline_path']} → {source['path']}"
                rows.append(
                    f"| <code>{markdown_escape( display_path )}</code> | "
                    f"{compact_metric_text( branch_change['baseline'] )} | "
                    f"{compact_metric_text( branch_change['candidate'] )} | "
                    f"{signed_integer( branch_change['not_covered_delta'] )} | "
                    f"{signed_integer( line_change['not_covered_delta'] )} |"
                )
            rows.extend(
                (
                    "",
                    "This is whole-file canonical context, ranked by newly missed "
                    "branches; the changed-code table below is the hunk-specific signal.",
                )
            )
    else:
        rows.extend(
            (
                "| Comparison | Result | Reason |",
                "|---|:---:|---|",
                "| Exact base → candidate | ℹ️ N/A | "
                f"{markdown_escape( aggregate['reason'] )} |",
            )
        )

    patch_metrics = patch["metrics"]
    rows.extend(
        (
            "",
            "#### Changed-code coverage",
            "",
            "| Signal | Covered | Missed | Total | Coverage | Interpretation |",
            "|---|---:|---:|---:|---:|---|",
            f"| Executable changed lines | "
            f"{patch_metrics['executable_lines']['covered']:,} | "
            f"{patch_metrics['executable_lines']['not_covered']:,} | "
            f"{patch_metrics['executable_lines']['total']:,} | "
            f"{percent_text( patch_metrics['executable_lines'] )} | "
            "Exact LLVM line mapping on added/modified head lines |",
            f"| Changed function start sites | "
            f"{patch_metrics['function_start_sites']['covered']:,} | "
            f"{patch_metrics['function_start_sites']['not_covered']:,} | "
            f"{patch_metrics['function_start_sites']['total']:,} | "
            f"{percent_text( patch_metrics['function_start_sites'] )} | "
            "Advisory; C++ instantiations may share a source definition |",
            f"| Native changed branch outcomes | "
            f"{patch_metrics['native_branch_outcomes']['covered']:,} | "
            f"{patch_metrics['native_branch_outcomes']['not_covered']:,} | "
            f"{patch_metrics['native_branch_outcomes']['total']:,} | "
            f"{percent_text( patch_metrics['native_branch_outcomes'] )} | "
            "Advisory and instantiation-weighted; not the canonical branch gate |",
            "",
            f"Git identified **{patch['changed_lines']:,}** added/modified head lines "
            f"across **{patch['changed_files']:,}** production files. "
            f"**{patch['non_instrumented_changed_lines']:,}** changed lines have no "
            "LLVM line mapping; these may be comments, declarations, braces, or "
            "configuration-dependent code.",
        )
    )
    if patch["excluded_changed_files"]:
        rows.append(
            f"The production coverage contract excludes an additional "
            f"**{patch['excluded_changed_lines']:,}** changed lines across "
            f"**{patch['excluded_changed_files']:,}** test/example/support files."
        )

    impacted_files = patch["top_impacted_files"]
    if impacted_files:
        rows.extend(
            (
                "",
                "#### Changed files with the largest coverage gaps",
                "",
                "| Source file | Change | Changed lines | Executable lines | "
                "Function starts | Native branches |",
                "|---|:---:|---:|---:|---:|---:|",
            )
        )
        for file_result in impacted_files:
            metrics = file_result["metrics"]
            rows.append(
                f"| <code>{markdown_escape( file_result['path'] )}</code> | "
                f"{file_result['status']} | {file_result['changed_lines']:,} | "
                f"{compact_metric_text( metrics['executable_lines'] )} | "
                f"{compact_metric_text( metrics['function_start_sites'] )} | "
                f"{compact_metric_text( metrics['native_branch_outcomes'] )} |"
            )

    impacted_locations = patch["top_impacted_locations"]
    rows.extend( ( "", "#### Highest-impact changed locations", "" ) )
    if impacted_locations:
        rows.extend(
            (
                "| Location | Coverage gap | Suggested test |",
                "|---|---|---|",
            )
        )
        for location in impacted_locations:
            gaps = []
            if location["executable_line"] and not location["line_covered"]:
                gaps.append( "executable line not hit" )
            missed_branches = location["native_branch_outcomes"]["not_covered"]
            if missed_branches:
                gaps.append(
                    f"{missed_branches:,} native branch outcome"
                    f"{'s' if missed_branches != 1 else ''} missed"
                )
            if location["function_start_covered"] is False:
                function_detail = "function start not entered"
                if location["function_names"]:
                    names = ", ".join(
                        f"<code>{markdown_escape( name )}</code>"
                        for name in location["function_names"]
                    )
                    function_detail += f" ({names}"
                    if location["additional_function_names"]:
                        function_detail += (
                            f", +{location['additional_function_names']:,} more"
                        )
                    function_detail += ")"
                gaps.append( function_detail )
            rows.append(
                f"| <code>{markdown_escape( location['path'] )}:"
                f"{location['line']}</code> | {'; '.join( gaps )} | "
                f"{markdown_escape( location['hint'] )} |"
            )
    else:
        rows.append(
            "No uncovered instrumented changed locations were found. If the patch "
            "has no executable changed lines, this section is informationally empty."
        )

    rows.extend(
        (
            "",
            "Changed-line coverage is exact. Function-start and native branch details "
            "are diagnostics because LLVM 20's detailed export retains C++ template "
            "and inline instantiations.",
            "",
            "This report is advisory and does not enforce a threshold.",
        )
    )
    return "\n".join( rows ) + "\n"


def analyze(
    repository: pathlib.Path,
    base_sha: str,
    head_sha: str,
    native_info: pathlib.Path,
    candidate_summary_path: pathlib.Path,
    baseline_summary_path: pathlib.Path | None = None,
) -> tuple[dict, str]:
    repository = validate_repository( repository )
    base_sha = validate_commit( repository, base_sha, "--base-sha" )
    head_sha = validate_commit( repository, head_sha, "--head-sha" )
    base_tree = commit_tree( repository, base_sha )
    head_tree = commit_tree( repository, head_sha )
    native_info = checked_input_file( native_info, "--native-info" )
    candidate_summary_path = checked_input_file(
        candidate_summary_path, "--candidate-summary"
    )
    candidate = validate_summary(
        read_json( candidate_summary_path, "candidate summary" ),
        "candidate summary",
    )
    if candidate["scope"] != PRODUCTION_SCOPE:
        raise ValueError(
            f"candidate summary scope {candidate['scope']!r} does not match "
            f"{PRODUCTION_SCOPE!r}"
        )

    baseline = None
    baseline_error = None
    if baseline_summary_path is not None:
        try:
            baseline_summary_path = checked_input_file(
                baseline_summary_path, "--baseline-summary"
            )
            baseline = validate_summary(
                read_json( baseline_summary_path, "baseline summary" ),
                "baseline summary",
            )
        except ValueError as error:
            baseline_error = f"N/A: baseline summary is invalid: {error}"

    changed_files = changed_line_ranges( repository, base_sha, head_sha )
    try:
        exclusion_pattern = re.compile( candidate["excluded_regex"] )
    except re.error as error:
        raise ValueError(
            f"candidate summary excluded_regex is invalid: {error}"
        ) from error
    included_changed_files = []
    excluded_changed_files = []
    for changed_file in changed_files:
        absolute_source = ( repository / changed_file["path"] ).as_posix()
        if exclusion_pattern.search( absolute_source ) is not None:
            excluded_changed_files.append( changed_file )
        else:
            included_changed_files.append( changed_file )
    changed_files = included_changed_files
    lcov = parse_native_lcov( native_info, repository, changed_files )
    patch = build_patch_report( changed_files, lcov )
    patch["excluded_changed_files"] = len( excluded_changed_files )
    patch["excluded_changed_lines"] = sum(
        len( changed_file["lines"] ) for changed_file in excluded_changed_files
    )
    aggregate = build_aggregate_comparison(
        candidate,
        baseline,
        baseline_error,
        base_sha,
        head_sha,
        base_tree,
        head_tree,
        changed_files,
    )
    report = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "scope": PRODUCTION_SCOPE,
        "base_sha": base_sha,
        "head_sha": head_sha,
        "base_tree_sha": base_tree,
        "head_tree_sha": head_tree,
        "candidate_summary_schema_version": candidate["schema_version"],
        "aggregate_comparison": aggregate,
        "patch": patch,
    }
    return report, build_markdown( report )


def checked_output_path(
    path: pathlib.Path,
    inputs: set[pathlib.Path],
    description: str,
) -> pathlib.Path:
    if path.is_symlink():
        raise ValueError( f"{description} must not be a symbolic link: {path}" )
    absolute = path.absolute()
    if absolute in inputs:
        raise ValueError( f"{description} must not overwrite an input file" )
    return absolute


def write_atomic( path: pathlib.Path, contents: str ) -> None:
    path.parent.mkdir( parents=True, exist_ok=True )
    if path.is_symlink():
        raise ValueError( f"refusing to replace symbolic link: {path}" )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary_path = pathlib.Path( temporary_name )
    try:
        # The CI analyzer runs as root in Docker while the GitHub runner reads
        # these files from a bind mount. mkstemp defaults to 0600, so make the
        # final compact report deliberately host-readable before replacement.
        os.fchmod( descriptor, 0o644 )
        with os.fdopen( descriptor, "w", encoding="utf-8" ) as stream:
            stream.write( contents )
        os.replace( temporary_path, path )
    except BaseException:
        temporary_path.unlink( missing_ok=True )
        raise


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report advisory source coverage for a Git pull-request diff."
    )
    parser.add_argument( "--repository", required=True, type=pathlib.Path )
    parser.add_argument( "--base-sha", required=True )
    parser.add_argument( "--head-sha", required=True )
    parser.add_argument( "--native-info", required=True, type=pathlib.Path )
    parser.add_argument( "--candidate-summary", required=True, type=pathlib.Path )
    parser.add_argument( "--baseline-summary", type=pathlib.Path )
    parser.add_argument( "--output-json", type=pathlib.Path )
    parser.add_argument( "--output-markdown", type=pathlib.Path )
    args = parser.parse_args()

    try:
        input_paths = {
            checked_input_file( args.native_info, "--native-info" ),
            checked_input_file( args.candidate_summary, "--candidate-summary" ),
        }
        if args.baseline_summary is not None:
            # A malformed baseline remains advisory, but an output must never be
            # allowed to overwrite the caller-provided path.
            input_paths.add( args.baseline_summary.absolute() )
        output_json = (
            checked_output_path( args.output_json, input_paths, "--output-json" )
            if args.output_json is not None
            else None
        )
        output_markdown = (
            checked_output_path(
                args.output_markdown, input_paths, "--output-markdown"
            )
            if args.output_markdown is not None
            else None
        )
        if output_json is not None and output_json == output_markdown:
            raise ValueError( "--output-json and --output-markdown must differ" )

        report, markdown = analyze(
            args.repository,
            args.base_sha,
            args.head_sha,
            args.native_info,
            args.candidate_summary,
            args.baseline_summary,
        )
        json_text = json.dumps(
            report, indent=2, sort_keys=True, ensure_ascii=False
        ) + "\n"
        if output_json is not None:
            write_atomic( output_json, json_text )
        if output_markdown is not None:
            write_atomic( output_markdown, markdown )
        print( markdown, end="" )
        return 0
    except ( OSError, ValueError ) as error:
        print( f"PR coverage comparison error: {error}", file=sys.stderr )
        return 2


if __name__ == "__main__":
    sys.exit( main() )
