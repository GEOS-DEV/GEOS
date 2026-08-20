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
CONSOLE_METRICS = ( "lines", "functions", "branches" )


def render_console_summary( summary: dict, thresholds: dict ) -> str:
    """Render a plain-text table for the GitHub Actions log."""
    rows = []
    all_passed = True
    for name in CONSOLE_METRICS:
        metric = summary["metrics"][name]
        threshold = thresholds[name]
        passed = metric["covered"] * 10000 >= threshold * metric["total"]
        all_passed = all_passed and passed
        rows.append(
            (
                DISPLAY_NAMES[name],
                f"{metric['percent']:.2f}%",
                f"{metric['covered']:,} / {metric['total']:,}",
                f"{metric['not_covered']:,}",
                f">= {threshold / 100:.2f}%",
                "✅ Pass" if passed else "❌ Fail",
            )
        )

    headers = ( "Metric", "Coverage", "Covered / total", "Missed", "Requirement", "Result" )
    widths = [
        max( len( headers[index] ), *( len( row[index] ) for row in rows ) )
        for index in range( len( headers ) )
    ]

    def border() -> str:
        return "+" + "+".join( "-" * ( width + 2 ) for width in widths ) + "+"

    def row( values ) -> str:
        return "|" + "|".join(
            f" {value:<{width}} "
            for value, width in zip( values, widths )
        ) + "|"

    status = "✅ PASS" if all_passed else "❌ FAIL"
    output = [
        "",
        "LLVM source coverage policy",
        f"Overall result: {status}",
        border(),
        row( headers ),
        border(),
        *( row( values ) for values in rows ),
        border(),
        "",
    ]
    return "\n".join( output )


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
    parser.add_argument( "--console", action="store_true" )
    parser.add_argument( "--tests-result" )
    parser.add_argument( "--tests-detail" )
    arguments = parser.parse_args()

    try:
        raw_summary = read_json( arguments.coverage_summary, "coverage summary" )
        raw_policy = read_json( arguments.thresholds, "coverage policy" )
        summary = coverage_contract.validate_summary( raw_summary )
        thresholds = coverage_contract.validate_policy( raw_policy, summary["scope"] )
        if arguments.console:
            print( render_console_summary( summary, thresholds ), end="" )
            return 0
        if arguments.tests_result is None or arguments.tests_detail is None:
            raise ValueError( "--tests-result and --tests-detail are required" )
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
