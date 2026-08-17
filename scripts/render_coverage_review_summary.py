#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

"""Render the concise reviewer-facing portion of the LLVM coverage summary."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


SCRIPT_DIRECTORY = pathlib.Path( __file__ ).resolve().parent
if str( SCRIPT_DIRECTORY ) not in sys.path:
    sys.path.insert( 0, str( SCRIPT_DIRECTORY ) )
import check_coverage_thresholds as coverage_contract


DISPLAY_NAMES = {
    "regions": "Regions",
    "functions": "Functions",
    "lines": "Lines",
    "branches": "Canonical branches",
}


def read_json( path: pathlib.Path, description: str ) -> dict:
    try:
        with path.open( encoding="utf-8" ) as stream:
            value = json.load( stream )
    except ( OSError, UnicodeError, json.JSONDecodeError ) as error:
        raise ValueError( f"cannot read {description} {path}: {error}" ) from error
    if not isinstance( value, dict ):
        raise ValueError( f"{description} must contain an object" )
    return value


def metric_text( metric: dict ) -> str:
    return (
        f"{metric['percent']:.2f}% "
        f"({metric['covered']:,}/{metric['total']:,})"
    )


def policy_row( summary: dict, thresholds: dict ) -> str:
    details = []
    all_passed = True
    for name in coverage_contract.DISPLAY_METRICS:
        if name not in thresholds:
            continue
        metric = summary["metrics"][name]
        threshold = thresholds[name]
        passed = metric["covered"] * 10000 >= threshold * metric["total"]
        all_passed = all_passed and passed
        result = "✅" if passed else "❌"
        details.append(
            f"{result} {DISPLAY_NAMES[name]} {metric_text( metric )} "
            f"≥ {threshold / 100:.2f}%"
        )
    status = "✅ Pass" if all_passed else "❌ Fail"
    return (
        f"| Enforced coverage policy | {status} | "
        f"{'<br>'.join( details )} |"
    )


def comparison_row( report: dict | None ) -> str:
    if report is None:
        return (
            "| Baseline comparison | ⏭️ Not applicable | "
            "This is not a pull-request coverage run. |"
        )

    aggregate = report.get( "aggregate_comparison", {} )
    if not aggregate.get( "comparable", False ):
        reason = aggregate.get( "reason", "No exact-base coverage summary was provided." )
        if reason == "No exact-base coverage summary was provided.":
            reason = "No exact-base baseline artifact was available."
        return f"| Baseline comparison | ℹ️ Unavailable | {reason} |"

    deltas = []
    for name in coverage_contract.DISPLAY_METRICS:
        comparison = aggregate["metrics"].get( name )
        if comparison is not None:
            deltas.append( comparison["percentage_point_delta"] )
    if all( delta >= 0 for delta in deltas ):
        result = "✅ Maintained/improved"
    elif all( delta <= 0 for delta in deltas ):
        result = "⚠️ Decreased"
    else:
        result = "⚠️ Mixed"
    details = "; ".join(
        f"{DISPLAY_NAMES[name]} {aggregate['metrics'][name]['percentage_point_delta']:+.2f} pp"
        for name in coverage_contract.DISPLAY_METRICS
        if name in aggregate["metrics"]
    )
    return f"| Baseline comparison | {result} | {details} |"


def changed_lines_row( report: dict | None ) -> str:
    if report is None:
        return (
            "| Changed production lines | ⏭️ Not applicable | "
            "This is not a pull-request coverage run. |"
        )

    metric = report["patch"]["metrics"]["executable_lines"]
    if metric["total"] == 0:
        return (
            "| Changed production lines | ℹ️ None | "
            "The production diff contains no executable changed lines. |"
        )

    locations = []
    for location in report["patch"].get( "top_impacted_locations", [] ):
        if location.get( "executable_line" ) and not location.get( "line_covered" ):
            locations.append( f"{location['path']}:{location['line']}" )
        if len( locations ) == 3:
            break
    details = (
        f"{metric['covered']:,}/{metric['total']:,} executable changed lines covered"
    )
    if locations:
        details += "; review " + ", ".join( f"<code>{path}</code>" for path in locations )
    if metric["not_covered"]:
        return f"| Changed production lines | ⚠️ Review | {details} |"
    return f"| Changed production lines | ✅ Covered | {details} |"


def render_summary(
    summary: dict,
    thresholds: dict,
    tests_result: str,
    tests_detail: str,
    report: dict | None,
) -> str:
    rows = [
        "### Reviewer summary",
        "",
        "| Signal | Result | Details |",
        "|---|:---:|---|",
        f"| Coverage smoke tests | {tests_result} | {tests_detail} |",
        policy_row( summary, thresholds ),
        comparison_row( report ),
        changed_lines_row( report ),
        "",
        "The policy row is the enforced gate. Baseline and changed-line signals "
        "are reviewer guidance; detailed diagnostics appear below.",
        "",
    ]
    return "\n".join( rows )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render the concise LLVM coverage review summary."
    )
    parser.add_argument( "--coverage-summary", required=True, type=pathlib.Path )
    parser.add_argument( "--thresholds", required=True, type=pathlib.Path )
    parser.add_argument( "--pr-report", type=pathlib.Path )
    parser.add_argument( "--tests-result", required=True )
    parser.add_argument( "--tests-detail", required=True )
    arguments = parser.parse_args()

    try:
        raw_summary = read_json( arguments.coverage_summary, "coverage summary" )
        raw_policy = read_json( arguments.thresholds, "coverage policy" )
        summary = coverage_contract.validate_summary( raw_summary )
        thresholds = coverage_contract.validate_policy( raw_policy, summary["scope"] )
        report = (
            read_json( arguments.pr_report, "pull-request coverage report" )
            if arguments.pr_report is not None
            else None
        )
        print(
            render_summary(
                summary,
                thresholds,
                arguments.tests_result,
                arguments.tests_detail,
                report,
            ),
            end="",
        )
        return 0
    except ( OSError, ValueError ) as error:
        print( f"Coverage review summary error: {error}", file=sys.stderr )
        return 1


if __name__ == "__main__":
    raise SystemExit( main() )
