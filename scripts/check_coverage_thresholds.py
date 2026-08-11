#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

"""Validate and enforce repository-owned LLVM coverage thresholds."""

import argparse
import html
import json
import math
import pathlib
import sys


COVERAGE_SUMMARY_SCHEMA_VERSION = 2
COVERAGE_POLICY_SCHEMA_VERSION = 1
CANONICAL_METRICS = ( "regions", "functions", "lines", "branches" )
DISPLAY_METRICS = ( "lines", "functions", "branches", "regions" )
SUPPLEMENTAL_METRICS = (
    ( "native_branch_outcomes", "Native branch outcomes" ),
)
if set( DISPLAY_METRICS ) != set( CANONICAL_METRICS ):
    raise RuntimeError( "display metrics must match canonical coverage metrics" )


def require_schema_version(
    document: dict, description: str, expected_version: int
) -> None:
    version = document.get( "schema_version" )
    if (
        isinstance( version, bool )
        or not isinstance( version, int )
        or version != expected_version
    ):
        raise ValueError(
            f"{description} schema_version must be the integer {expected_version}"
        )


def load_json( path: pathlib.Path ) -> dict:
    try:
        with path.open( encoding="utf-8" ) as stream:
            document = json.load( stream )
    except ( OSError, json.JSONDecodeError ) as error:
        raise ValueError( f"cannot read {path}: {error}" ) from error
    if not isinstance( document, dict ):
        raise ValueError( f"{path} must contain a JSON object" )
    return document


def require_integer( value, description: str, minimum: int = 0 ) -> int:
    if isinstance( value, bool ) or not isinstance( value, int ) or value < minimum:
        raise ValueError( f"{description} must be an integer >= {minimum}" )
    return value


def validate_metric( metric: object, description: str ) -> dict:
    if not isinstance( metric, dict ):
        raise ValueError( f"{description} must be an object" )

    covered = require_integer( metric.get( "covered" ), f"{description}.covered" )
    total = require_integer( metric.get( "total" ), f"{description}.total", 1 )
    not_covered = require_integer(
        metric.get( "not_covered" ), f"{description}.not_covered"
    )
    if covered > total:
        raise ValueError( f"{description}.covered cannot exceed total" )
    if not_covered != total - covered:
        raise ValueError( f"{description}.not_covered is inconsistent" )

    percent = metric.get( "percent" )
    if isinstance( percent, bool ) or not isinstance( percent, ( int, float ) ):
        raise ValueError( f"{description}.percent must be numeric" )
    if not math.isfinite( percent ):
        raise ValueError( f"{description}.percent must be finite" )
    expected_percent = round( 100.0 * covered / total, 6 )
    if abs( float( percent ) - expected_percent ) > 1.0e-6:
        raise ValueError( f"{description}.percent is inconsistent with its counts" )

    return {
        "covered": covered,
        "total": total,
        "not_covered": not_covered,
        "percent": expected_percent,
    }


def validate_branch_gaps( raw_gaps: object, scope: str ) -> list[dict]:
    if not isinstance( raw_gaps, list ) or len( raw_gaps ) > 5:
        raise ValueError( "top_branch_gaps must be an array with at most 5 entries" )

    scope_prefix = f"{scope.rstrip( '/' )}/"
    gaps = []
    seen_paths = set()
    for index, raw_gap in enumerate( raw_gaps ):
        if not isinstance( raw_gap, dict ):
            raise ValueError( f"top_branch_gaps[{index}] must be an object" )
        path = raw_gap.get( "path" )
        if not isinstance( path, str ) or not path.startswith( scope_prefix ):
            raise ValueError( f"top_branch_gaps[{index}].path is outside the scope" )
        normalized_path = pathlib.PurePosixPath( path )
        if (
            len( path ) == len( scope_prefix )
            or normalized_path.as_posix() != path
            or ".." in normalized_path.parts
            or path in seen_paths
        ):
            raise ValueError( "top_branch_gaps paths must be normalized and unique" )
        seen_paths.add( path )

        gap = validate_metric( raw_gap, f"top_branch_gaps[{index}]" )
        if gap["not_covered"] == 0:
            raise ValueError( "top_branch_gaps entries must have missed branches" )
        gaps.append( { "path": path, **gap } )

    expected_order = sorted(
        gaps, key=lambda gap: ( -gap["not_covered"], gap["path"] )
    )
    if gaps != expected_order:
        raise ValueError( "top_branch_gaps must be sorted by missed branches" )
    return gaps


def validate_summary( document: dict ) -> dict:
    require_schema_version(
        document, "coverage summary", COVERAGE_SUMMARY_SCHEMA_VERSION
    )
    scope = document.get( "scope" )
    if not isinstance( scope, str ) or not scope:
        raise ValueError( "coverage summary scope must be a nonempty string" )

    excluded_regex = document.get( "excluded_regex" )
    if not isinstance( excluded_regex, str ) or not excluded_regex:
        raise ValueError( "coverage summary excluded_regex must be nonempty" )

    tool = document.get( "tool" )
    if not isinstance( tool, dict ) or tool.get( "name" ) != "llvm-cov":
        raise ValueError( "coverage summary tool.name must be llvm-cov" )
    tool_major = require_integer( tool.get( "major" ), "tool.major", 1 )

    raw_inputs = document.get( "inputs" )
    if not isinstance( raw_inputs, dict ):
        raise ValueError( "coverage summary inputs must be an object" )
    inputs = {
        "profiles": require_integer( raw_inputs.get( "profiles" ), "inputs.profiles", 1 ),
        "coverage_objects": require_integer(
            raw_inputs.get( "coverage_objects" ), "inputs.coverage_objects", 1
        ),
        "zero_hash_mappings": require_integer(
            raw_inputs.get( "zero_hash_mappings" ), "inputs.zero_hash_mappings"
        ),
    }

    raw_metrics = document.get( "metrics" )
    if not isinstance( raw_metrics, dict ):
        raise ValueError( "coverage summary metrics must be an object" )
    metrics = {
        name: validate_metric( raw_metrics.get( name ), f"metrics.{name}" )
        for name in CANONICAL_METRICS
    }

    raw_supplemental = document.get( "supplemental" )
    if not isinstance( raw_supplemental, dict ):
        raise ValueError( "coverage summary supplemental must be an object" )
    supplemental = {
        name: validate_metric( raw_supplemental.get( name ), f"supplemental.{name}" )
        for name, _ in SUPPLEMENTAL_METRICS
    }
    branch_gaps = validate_branch_gaps( document.get( "top_branch_gaps" ), scope )
    if sum( gap["total"] for gap in branch_gaps ) > metrics["branches"]["total"]:
        raise ValueError( "top_branch_gaps totals exceed canonical branch coverage" )
    return {
        "scope": scope,
        "excluded_regex": excluded_regex,
        "tool": { "name": "llvm-cov", "major": tool_major },
        "inputs": inputs,
        "metrics": metrics,
        "top_branch_gaps": branch_gaps,
        "supplemental": supplemental,
    }


def validate_policy( document: dict, expected_scope: str ) -> dict:
    require_schema_version(
        document, "coverage policy", COVERAGE_POLICY_SCHEMA_VERSION
    )
    if document.get( "scope" ) != expected_scope:
        raise ValueError(
            f"coverage policy scope {document.get('scope')!r} does not match "
            f"summary scope {expected_scope!r}"
        )

    raw_thresholds = document.get( "minimum_basis_points" )
    if not isinstance( raw_thresholds, dict ) or not raw_thresholds:
        raise ValueError( "coverage policy minimum_basis_points must be nonempty" )

    thresholds = {}
    for name, value in raw_thresholds.items():
        if name not in CANONICAL_METRICS:
            raise ValueError( f"unsupported coverage metric in policy: {name}" )
        thresholds[name] = require_integer(
            value, f"minimum_basis_points.{name}"
        )
        if thresholds[name] > 10000:
            raise ValueError( f"minimum_basis_points.{name} cannot exceed 10000" )
    return thresholds


def percent_text( covered: int, total: int ) -> str:
    return f"{100.0 * covered / total:.2f}%"


def build_markdown( summary: dict, thresholds: dict ) -> tuple[str, bool]:
    threshold_results = {
        name: summary["metrics"][name]["covered"] * 10000
        >= threshold * summary["metrics"][name]["total"]
        for name, threshold in thresholds.items()
    }
    all_passed = all( threshold_results.values() )
    display_names = {
        "regions": "Regions",
        "functions": "Functions",
        "lines": "Lines",
        "branches": "Canonical branches",
    }
    metric_rows = []
    for name in DISPLAY_METRICS:
        metric = summary["metrics"][name]
        threshold = thresholds.get( name )
        if threshold is None:
            requirement = "—"
            result = "ℹ️ Measured"
        else:
            passed = threshold_results[name]
            requirement = f"≥ {threshold / 100:.2f}%"
            margin = 100.0 * metric["covered"] / metric["total"] - threshold / 100
            if 0 < abs( margin ) < 0.01:
                direction = "above" if margin > 0 else "below"
                margin_display = f"<0.01 pp {direction}"
            else:
                margin_display = f"{margin:+.2f} pp"
            result = f"{'✅ Pass' if passed else '❌ Fail'} ({margin_display})"
        metric_rows.append(
            f"| {display_names[name]} | {metric['covered']:,} | "
            f"{metric['not_covered']:,} | {metric['total']:,} | "
            f"{percent_text(metric['covered'], metric['total'])} | "
            f"{requirement} | {result} |"
        )

    policy_result = "passed" if all_passed else "failed"
    policy_icon = "✅" if all_passed else "❌"
    rows = [
        f"### {policy_icon} Coverage policy {policy_result}",
        "",
        f"Production scope: `{summary['scope']}` · "
        f"Reporter: `{summary['tool']['name']} {summary['tool']['major']}`",
        "",
        "| Metric | Covered | Missed | Total | Coverage | Requirement | Result |",
        "|---|---:|---:|---:|---:|---:|:---:|",
        *metric_rows,
        "",
        "Percentages are rounded for display; policy thresholds use exact covered/total counts.",
    ]
    branch_gaps = summary["top_branch_gaps"]
    if branch_gaps:
        rows.extend(
            (
                "",
                "### Largest branch gaps",
                "",
                "| Source file | Covered | Missed | Total | Coverage |",
                "|---|---:|---:|---:|---:|",
            )
        )
        for gap in branch_gaps:
            scope_prefix = f"{summary['scope'].rstrip( '/' )}/"
            display_path = html.escape( gap["path"][len( scope_prefix ):] )
            display_path = display_path.replace( "|", "&#124;" )
            display_path = display_path.replace( "\n", "&#10;" ).replace(
                "\r", "&#13;"
            )
            rows.append(
                f"| <code>{display_path}</code> | {gap['covered']:,} | "
                f"{gap['not_covered']:,} | {gap['total']:,} | "
                f"{percent_text(gap['covered'], gap['total'])} |"
            )
        rows.extend(
            (
                "",
                "Ranked by missed canonical branch outcomes within the production scope.",
            )
        )

    rows.extend(
        (
            "",
            "### Coverage collection",
            "",
            "| Input | Value |",
            "|---|---:|",
            f"| Profile files merged | {summary['inputs']['profiles']:,} |",
            f"| Instrumented product objects | {summary['inputs']['coverage_objects']:,} |",
            f"| Audited zero-hash mappings | {summary['inputs']['zero_hash_mappings']:,} |",
            "",
            "Zero-hash mappings have no profile counters and are audited during report "
            "generation; any nonzero or unaudited mismatch fails the job.",
            "",
            "<details>",
            "<summary>Branch diagnostic (not gated)</summary>",
            "",
            "| Diagnostic | Covered | Total | Coverage |",
            "|---|---:|---:|---:|",
        )
    )
    for name, display_name in SUPPLEMENTAL_METRICS:
        metric = summary["supplemental"][name]
        rows.append(
            f"| {display_name} | {metric['covered']:,} | "
            f"{metric['total']:,} | {percent_text(metric['covered'], metric['total'])} |"
        )

    rows.extend(
        (
            "",
            "Native outcomes are instantiation-weighted: C++ template and inline-function "
            "instances may appear more than once. This diagnostic is not the policy metric.",
            "",
            "</details>",
        )
    )
    return "\n".join( rows ) + "\n", all_passed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument( "summary", type=pathlib.Path )
    parser.add_argument( "policy", type=pathlib.Path )
    parser.add_argument( "--markdown", type=pathlib.Path )
    args = parser.parse_args()

    if args.markdown is not None:
        try:
            args.markdown.unlink( missing_ok=True )
        except OSError as error:
            print(
                f"coverage threshold error: cannot clear {args.markdown}: {error}",
                file=sys.stderr,
            )
            return 2

    try:
        summary = validate_summary( load_json( args.summary ) )
        thresholds = validate_policy( load_json( args.policy ), summary["scope"] )
        markdown, passed = build_markdown( summary, thresholds )
    except ValueError as error:
        print( f"coverage threshold error: {error}", file=sys.stderr )
        return 2

    print( markdown, end="" )
    if args.markdown is not None:
        args.markdown.parent.mkdir( parents=True, exist_ok=True )
        args.markdown.write_text( markdown, encoding="utf-8" )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit( main() )
