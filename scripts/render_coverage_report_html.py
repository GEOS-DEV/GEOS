#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

"""Render the LLVM source-coverage artifact as one self-contained HTML page."""

from __future__ import annotations

import argparse
import html
import json
import os
import pathlib
import re
import sys
import tempfile


SCRIPT_DIRECTORY = pathlib.Path( __file__ ).resolve().parent
if str( SCRIPT_DIRECTORY ) not in sys.path:
    sys.path.insert( 0, str( SCRIPT_DIRECTORY ) )
import check_coverage_thresholds as coverage_contract


METRIC_NAMES = ( "lines", "functions", "branches", "regions" )
METRIC_LABELS = {
    "branches": "Canonical branches",
    "functions": "Functions",
    "lines": "Lines",
    "regions": "Regions",
}
PATCH_METRIC_LABELS = {
    "executable_lines": "Executable changed lines",
    "function_start_sites": "Changed function start sites",
    "native_branch_outcomes": "Native changed branch outcomes",
}
DETAIL_FILES = (
    ( "CTest output", "ctest.log" ),
    ( "LLVM mapping-integrity log", "mapping-integrity.log" ),
    ( "Native export log", "native-export.log" ),
    ( "Summary export log", "summary-export.log" ),
    ( "Branch summary", "branch-summary.txt" ),
    ( "Coverage objects", "coverage-objects.txt" ),
    ( "Toolchain", "toolchain.txt" ),
    ( "CMake configuration", "CMakeCache.txt" ),
    ( "Coverage status source", "coverage-status.md" ),
    ( "Coverage policy source", "coverage-summary.md" ),
    ( "Pull-request coverage source", "pr-coverage.md" ),
    ( "Baseline status", "coverage-baseline-status.json" ),
    ( "Coverage summary JSON", "coverage-summary.json" ),
    ( "Pull-request coverage JSON", "pr-coverage.json" ),
    ( "LLVM summary JSON", "llvm-summary.json" ),
)


def escape( value: object ) -> str:
    return html.escape( str( value ), quote=True )


def read_text( path: pathlib.Path ) -> str | None:
    if not path.is_file() or path.is_symlink():
        return None
    try:
        return path.read_text( encoding="utf-8", errors="replace" )
    except OSError:
        return None


def read_json( path: pathlib.Path ) -> dict | None:
    text = read_text( path )
    if text is None:
        return None
    try:
        document = json.loads( text )
    except json.JSONDecodeError:
        return None
    return document if isinstance( document, dict ) else None


def write_atomic( path: pathlib.Path, contents: str ) -> None:
    """Replace a report even when a container created the old file as root."""
    path.parent.mkdir( parents=True, exist_ok=True )
    if path.is_symlink():
        raise ValueError( f"refusing to replace symbolic link: {path}" )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary_path = pathlib.Path( temporary_name )
    try:
        os.fchmod( descriptor, 0o644 )
        with os.fdopen( descriptor, "w", encoding="utf-8" ) as stream:
            stream.write( contents )
            stream.flush()
            os.fsync( stream.fileno() )
        os.replace( temporary_path, path )
        temporary_path = None
    except BaseException:
        temporary_path.unlink( missing_ok=True )
        raise


def validated_summary( document: dict | None ) -> dict | None:
    if document is None:
        return None
    try:
        return coverage_contract.validate_summary( document )
    except ValueError:
        return document


def validated_policy( document: dict | None, scope: str | None ) -> dict:
    if document is None or scope is None:
        return {}
    try:
        return coverage_contract.validate_policy( document, scope )
    except ValueError:
        raw_thresholds = document.get( "minimum_basis_points", {} )
        if not isinstance( raw_thresholds, dict ):
            return {}
        return {
            name: value
            for name, value in raw_thresholds.items()
            if name in METRIC_NAMES and isinstance( value, int )
        }


def number( value: object ) -> str:
    if isinstance( value, int ) and not isinstance( value, bool ):
        return f"{value:,}"
    return escape( value )


def metric_percent( metric: dict | None ) -> str:
    if not isinstance( metric, dict ):
        return "N/A"
    total = metric.get( "total", 0 )
    covered = metric.get( "covered", 0 )
    if not isinstance( total, int ) or total == 0:
        return "N/A"
    return f"{100.0 * covered / total:.2f}%"


def metric_counts( metric: dict | None ) -> str:
    if not isinstance( metric, dict ):
        return "N/A"
    return f"{number( metric.get( 'covered', 0 ) )}/{number( metric.get( 'total', 0 ) )}"


def policy_result( metric: dict | None, threshold: int | None ) -> tuple[str, str]:
    if not isinstance( metric, dict ) or threshold is None:
        return "info", "Measured"
    total = metric.get( "total", 0 )
    covered = metric.get( "covered", 0 )
    if not isinstance( total, int ) or not isinstance( covered, int ) or total == 0:
        return "info", "N/A"
    passed = covered * 10000 >= threshold * total
    return ( "pass", "Pass" ) if passed else ( "fail", "Fail" )


def badge( state: str, label: str ) -> str:
    icons = { "pass": "✓", "fail": "✕", "warn": "!", "info": "i" }
    return (
        f'<span class="badge {escape( state )}">'
        f'<span class="badge-icon">{icons.get( state, "i" )}</span>'
        f"{escape( label )}</span>"
    )


def table( headers: list[str], rows: list[list[str]], class_name="" ) -> str:
    header_html = "".join( f"<th>{header}</th>" for header in headers )
    body = []
    for row in rows:
        body.append( "<tr>" + "".join( f"<td>{cell}</td>" for cell in row ) + "</tr>" )
    return (
        f'<div class="table-wrap"><table class="{escape( class_name )}">'
        f"<thead><tr>{header_html}</tr></thead><tbody>{''.join( body )}</tbody>"
        "</table></div>"
    )


def card( title: str, result: str, detail: str, state: str ) -> str:
    return (
        f'<section class="card {escape( state )}">'
        f"<div class=\"card-title\">{escape( title )}</div>"
        f'<div class="card-result">{badge( state, result )}</div>'
        f'<div class="card-detail">{detail}</div>'
        "</section>"
    )


def markdown_table_value( document: str, label: str ) -> tuple[str, str] | None:
    for line in document.splitlines():
        cells = [ cell.strip() for cell in line.strip().strip( "|" ).split( "|" ) ]
        if len( cells ) >= 3 and cells[0] == label:
            return cells[1], cells[2]
    return None


def status_details( status_markdown: str | None ) -> dict:
    status = "info"
    result = "Unavailable"
    if status_markdown:
        if re.search( r"^### ✅ LLVM source coverage passed", status_markdown, re.MULTILINE ):
            status, result = "pass", "Passed"
        elif re.search( r"^### ❌ LLVM source coverage failed", status_markdown, re.MULTILINE ):
            status, result = "fail", "Failed"
    tests = markdown_table_value( status_markdown or "", "Coverage smoke suite" )
    if tests is None:
        tests = markdown_table_value( status_markdown or "", "Coverage smoke CTest" )
    configuration = {}
    for label in ( "Runner", "Container", "Immutable image ID (short)", "Toolchain", "Build", "Test scheduling" ):
        value = markdown_table_value( status_markdown or "", label )
        if value is not None:
            configuration[label] = value[1].strip( "`" )
    return {
        "state": status,
        "result": result,
        "tests": tests,
        "configuration": configuration,
    }


def metric_rows( summary: dict | None, thresholds: dict ) -> list[list[str]]:
    if summary is None:
        return []
    rows = []
    metrics = summary.get( "metrics", {} )
    for name in METRIC_NAMES:
        metric = metrics.get( name )
        threshold = thresholds.get( name )
        state, result = policy_result( metric, threshold )
        requirement = "—" if threshold is None else f"≥ {threshold / 100:.2f}%"
        percent = metric_percent( metric )
        width = 0.0
        if isinstance( metric, dict ) and isinstance( metric.get( "percent" ), ( int, float ) ):
            width = max( 0.0, min( 100.0, float( metric["percent"] ) ) )
        progress = (
            f'<div class="progress"><span style="width: {width:.2f}%"></span></div>'
        )
        rows.append(
            [
                escape( METRIC_LABELS[name] ),
                f"<strong>{escape( percent )}</strong>{progress}",
                number( metric.get( "covered", 0 ) ) if isinstance( metric, dict ) else "N/A",
                number( metric.get( "not_covered", 0 ) ) if isinstance( metric, dict ) else "N/A",
                number( metric.get( "total", 0 ) ) if isinstance( metric, dict ) else "N/A",
                escape( requirement ),
                badge( state, result ),
            ]
        )
    return rows


def policy_card( summary: dict | None, thresholds: dict ) -> tuple[str, str, str]:
    if summary is None:
        return "info", "Unavailable", "No validated coverage summary was produced."
    results = [ policy_result( summary.get( "metrics", {} ).get( name ), threshold )[0]
                for name, threshold in thresholds.items() ]
    if not results:
        return "info", "Measured", "No repository thresholds are configured."
    if "fail" in results:
        return "fail", "Failed", "At least one enforced coverage floor is below the required value."
    if "info" in results:
        return "info", "Unavailable", "One or more configured coverage metrics are unavailable."
    return "pass", "Passed", "All configured lines, functions, and branch thresholds are satisfied."


def baseline_card( report: dict | None ) -> tuple[str, str, str]:
    if report is None:
        return "info", "Not applicable", "This artifact was not produced for a pull request."
    aggregate = report.get( "aggregate_comparison", {} )
    if not aggregate.get( "comparable", False ):
        return "info", "Unavailable", aggregate.get(
            "reason", "No exact-base coverage summary was available."
        )
    deltas = [
        comparison.get( "percentage_point_delta", 0.0 )
        for comparison in aggregate.get( "metrics", {} ).values()
    ]
    if deltas and all( delta >= 0 for delta in deltas ):
        return "pass", "Maintained / improved", "Every comparable repository metric stayed level or increased."
    if deltas and all( delta <= 0 for delta in deltas ):
        return "warn", "Decreased", "At least one comparable repository metric decreased."
    if not deltas:
        return "info", "Unavailable", "No comparable repository metrics were provided."
    return "warn", "Mixed", "Some comparable repository metrics increased while others decreased."


def changed_lines_card( report: dict | None ) -> tuple[str, str, str]:
    if report is None:
        return "info", "Not applicable", "This artifact was not produced for a pull request."
    metric = report.get( "patch", {} ).get( "metrics", {} ).get( "executable_lines", {} )
    total = metric.get( "total", 0 )
    missed = metric.get( "not_covered", 0 )
    if total == 0:
        return "info", "No executable changes", "The production diff contains no executable changed lines."
    if missed:
        return "warn", "Review needed", f"{number( metric.get( 'covered', 0 ) )}/{number( total )} executable changed lines are covered."
    return "pass", "Covered", f"All {number( total )} executable changed lines are covered."


def coverage_section( summary: dict | None, thresholds: dict ) -> str:
    if summary is None:
        return '<section id="coverage"><h2>Repository coverage</h2><p class="empty">Coverage metrics are unavailable.</p></section>'
    rows = metric_rows( summary, thresholds )
    scope = escape( summary.get( "scope", "unknown" ) )
    reporter = f"{summary.get( 'tool', {} ).get( 'name', 'unknown' )} {summary.get( 'tool', {} ).get( 'major', '?' )}"
    result = (
        f'<section id="coverage"><h2>Repository coverage</h2>'
        f'<p class="section-intro">Production scope <code>{scope}</code> · Reporter <code>{escape( reporter )}</code>. '
        "Thresholds are evaluated using exact covered and total counts.</p>"
        + table(
            [ "Metric", "Coverage", "Covered", "Missed", "Total", "Requirement", "Result" ],
            rows,
            "metrics",
        )
    )
    gaps = summary.get( "top_branch_gaps", [] )
    if gaps:
        gap_rows = []
        scope_prefix = summary.get( "scope", "" ).rstrip( "/" ) + "/"
        for gap in gaps:
            path = gap.get( "path", "" )
            if path.startswith( scope_prefix ):
                path = path[len( scope_prefix ):]
            gap_rows.append(
                [
                    f"<code>{escape( path )}</code>",
                    number( gap.get( "covered", 0 ) ),
                    number( gap.get( "not_covered", 0 ) ),
                    number( gap.get( "total", 0 ) ),
                    escape( metric_percent( gap ) ),
                ]
            )
        result += (
            "<h3>Largest branch gaps</h3>"
            + table( [ "Source file", "Covered", "Missed", "Total", "Coverage" ], gap_rows )
            + '<p class="caption">Ranked by missed canonical branch outcomes within the production scope.</p>'
        )
    supplemental = summary.get( "supplemental", {} ).get( "native_branch_outcomes" )
    inputs = summary.get( "inputs", {} )
    collection_rows = [
        [ "Profile files merged", number( inputs.get( "profiles", "N/A" ) ) ],
        [ "Instrumented product objects", number( inputs.get( "coverage_objects", "N/A" ) ) ],
        [ "Audited zero-hash mappings", number( inputs.get( "zero_hash_mappings", "N/A" ) ) ],
    ]
    result += (
        "<h3>Coverage collection</h3>"
        + table( [ "Input", "Value" ], collection_rows )
        + '<p class="caption">Zero-hash mappings have no profile counters and are audited during report generation; any nonzero or unaudited mismatch fails the job.</p>'
    )
    if supplemental is not None:
        result += (
            '<div class="diagnostic-callout"><strong>Diagnostic only:</strong> native branch outcomes '
            f"{escape( metric_counts( supplemental ) )} ({escape( metric_percent( supplemental ))}). "
            "They are instantiation-weighted and are not the canonical branch gate.</div>"
        )
    return result + "</section>"


def comparison_table( report: dict ) -> str:
    aggregate = report.get( "aggregate_comparison", {} )
    if not aggregate.get( "comparable", False ):
        reason = escape( aggregate.get( "reason", "No exact-base coverage summary was available." ) )
        return table( [ "Comparison", "Result", "Reason" ], [ [ "Exact base → candidate", badge( "info", "Unavailable" ), reason ] ] )
    rows = []
    for name in METRIC_NAMES:
        comparison = aggregate.get( "metrics", {} ).get( name )
        if comparison is None:
            continue
        delta = comparison.get( "percentage_point_delta", 0.0 )
        state = "pass" if delta >= 0 else "warn"
        rows.append(
            [
                escape( METRIC_LABELS[name] ),
                escape( metric_counts( comparison.get( "baseline" ) ) + " (" + metric_percent( comparison.get( "baseline" ) ) + ")" ),
                escape( metric_counts( comparison.get( "candidate" ) ) + " (" + metric_percent( comparison.get( "candidate" ) ) + ")" ),
                f'<span class="delta {state}">{escape( f"{delta:+.2f} pp" )}</span>',
                escape( f"{comparison.get( 'covered_delta', 0):+,}" ),
                escape( f"{comparison.get( 'not_covered_delta', 0):+,}" ),
                escape( f"{comparison.get( 'total_delta', 0):+,}" ),
            ]
        )
    return table(
        [ "Metric", "Exact base", "Candidate", "Coverage change", "Covered Δ", "Missed Δ", "Total Δ" ],
        rows,
    )


def pull_request_section( report: dict | None ) -> str:
    if report is None:
        return '<section id="pull-request"><h2>Pull-request coverage</h2><p class="empty">No pull-request coverage report was produced.</p></section>'
    aggregate = report.get( "aggregate_comparison", {} )
    patch = report.get( "patch", {} )
    base = escape( str( report.get( "base_sha", "unknown" ) )[:12] )
    head = escape( str( report.get( "head_sha", "unknown" ) )[:12] )
    scope = escape( report.get( "scope", "unknown" ) )
    result = (
        f'<section id="pull-request"><h2>Pull-request coverage</h2>'
        f'<p class="section-intro">Effective source diff <code>{base}</code> → <code>{head}</code> · Production scope <code>{scope}</code>.</p>'
        "<h3>Project coverage change</h3>"
        + comparison_table( report )
    )
    if aggregate.get( "comparable", False ):
        result += '<p class="caption">Comparison is exact because commit, tree, and measurement-contract provenance match.</p>'
    else:
        result += '<p class="caption">Baseline comparison is advisory and is unavailable when no exact-base artifact is present.</p>'

    patch_metrics = patch.get( "metrics", {} )
    patch_rows = []
    interpretations = {
        "executable_lines": "Exact LLVM line mapping on added/modified head lines",
        "function_start_sites": "Advisory; C++ instantiations may share a source definition",
        "native_branch_outcomes": "Advisory and instantiation-weighted; not the canonical branch gate",
    }
    for name, label in PATCH_METRIC_LABELS.items():
        metric = patch_metrics.get( name, {} )
        state = "pass" if metric.get( "total", 0 ) == 0 or metric.get( "not_covered", 0 ) == 0 else "warn"
        patch_rows.append(
            [
                escape( label ),
                number( metric.get( "covered", 0 ) ),
                number( metric.get( "not_covered", 0 ) ),
                number( metric.get( "total", 0 ) ),
                escape( metric_percent( metric ) ),
                badge( state, "Covered" if state == "pass" else "Review" ),
                escape( interpretations[name] ),
            ]
        )
    result += "<h3>Changed-code coverage</h3>" + table(
        [ "Signal", "Covered", "Missed", "Total", "Coverage", "Review", "Interpretation" ],
        patch_rows,
    )
    result += (
        f'<p class="section-intro">Git identified <strong>{number( patch.get( "changed_lines", 0 ) )}</strong> '
        f"added/modified head lines across <strong>{number( patch.get( 'changed_files', 0 ) )}</strong> production files. "
        f"<strong>{number( patch.get( 'non_instrumented_changed_lines', 0 ) )}</strong> changed lines have no LLVM line mapping."
        "</p>"
    )

    impacted_files = patch.get( "top_impacted_files", [] )
    if impacted_files:
        file_rows = []
        for source in impacted_files:
            metrics = source.get( "metrics", {} )
            file_rows.append(
                [
                    f'<code>{escape( source.get( "path", "" ) )}</code>',
                    escape( source.get( "status", "" ) ),
                    number( source.get( "changed_lines", 0 ) ),
                    escape( metric_counts( metrics.get( "executable_lines" ) ) ),
                    escape( metric_counts( metrics.get( "function_start_sites" ) ) ),
                    escape( metric_counts( metrics.get( "native_branch_outcomes" ) ) ),
                ]
            )
        result += "<h3>Changed files with the largest coverage gaps</h3>" + table(
            [ "Source file", "Change", "Changed lines", "Executable lines", "Function starts", "Native branches" ],
            file_rows,
        )

    locations = patch.get( "top_impacted_locations", [] )
    result += "<h3>Highest-impact changed locations</h3>"
    if locations:
        location_rows = []
        for location in locations:
            gaps = []
            if location.get( "executable_line" ) and not location.get( "line_covered" ):
                gaps.append( "executable line not hit" )
            missed_branches = location.get( "native_branch_outcomes", {} ).get( "not_covered", 0 )
            if missed_branches:
                gaps.append( f"{number( missed_branches )} native branch outcome(s) missed" )
            if location.get( "function_start_covered" ) is False:
                gaps.append( "function start not entered" )
            location_rows.append(
                [
                    f'<code>{escape( location.get( "path", "" ) )}:{number( location.get( "line", "" ) )}</code>',
                    escape( "; ".join( gaps ) ),
                    escape( location.get( "hint", "" ) ),
                ]
            )
        result += table( [ "Location", "Coverage gap", "Suggested test" ], location_rows )
    else:
        result += '<p class="empty">No uncovered instrumented changed locations were found.</p>'
    return result + (
        '<p class="caption">Changed-line coverage is exact. Function-start and native branch details are diagnostics because LLVM 20 retains C++ template and inline-function instantiations.</p></section>'
    )


def details_section( artifact_dir: pathlib.Path ) -> str:
    sections = []
    for label, filename in DETAIL_FILES:
        path = artifact_dir / filename
        text = read_text( path )
        if text is None:
            continue
        if filename.endswith( ".json" ):
            document = read_json( path )
            if document is not None:
                text = json.dumps( document, indent=2, sort_keys=True, ensure_ascii=False ) + "\n"
        sections.append(
            f'<details><summary>{escape( label )} <span class="file-name">{escape( filename )}</span> '
            f'<span class="file-size">{number( len( text.encode( "utf-8" ) ) )} bytes</span></summary>'
            f'<pre>{escape( text )}</pre></details>'
        )
    if not sections:
        return '<section id="details"><h2>Diagnostics and raw data</h2><p class="empty">No diagnostic files were produced.</p></section>'
    return '<section id="details"><h2>Diagnostics and raw data</h2><p class="section-intro">All intermediate reports are retained below in collapsible sections, so this single HTML file remains the complete downloadable artifact.</p>' + "".join( sections ) + "</section>"


def render_html( artifact_dir: pathlib.Path ) -> str:
    status_markdown = read_text( artifact_dir / "coverage-status.md" )
    status = status_details( status_markdown )
    summary = validated_summary( read_json( artifact_dir / "coverage-summary.json" ) )
    policy_document = read_json( pathlib.Path( __file__ ).parents[1] / ".github" / "coverage-thresholds.json" )
    thresholds = validated_policy( policy_document, summary.get( "scope" ) if summary else None )
    report = read_json( artifact_dir / "pr-coverage.json" )
    policy_state, policy_result_text, policy_detail = policy_card( summary, thresholds )
    baseline_state, baseline_result, baseline_detail = baseline_card( report )
    changed_state, changed_result, changed_detail = changed_lines_card( report )
    test_result, test_detail = status.get( "tests" ) or ( "Unavailable", "No CTest summary was recorded." )
    overall_state = status["state"]
    overall_result = status["result"]
    if overall_state == "info":
        overall_state = "fail" if policy_state == "fail" else policy_state
        overall_result = "Failed" if overall_state == "fail" else "Available"
    summary_metrics = summary or {}
    scope = escape( summary_metrics.get( "scope", "src/coreComponents" ) )
    measurement = summary_metrics.get( "measurement", {} )
    commit = escape( str( measurement.get( "commit_sha", "unknown" ) )[:12] )
    configuration = status.get( "configuration", {} )
    config_rows = [ [ escape( label ), escape( value ) ] for label, value in configuration.items() ]
    if not config_rows and measurement:
        toolchain = measurement.get( "toolchain", {} )
        config_rows = [
            [ "Compiler target", escape( toolchain.get( "compiler_target", "unknown" ) ) ],
            [ "LLVM", escape( summary_metrics.get( "tool", {} ).get( "major", "unknown" ) ) ],
        ]
    nav = ''.join(
        f'<a href="#{anchor}">{escape( label )}</a>'
        for anchor, label in (
            ( "coverage", "Coverage" ),
            ( "pull-request", "Pull request" ),
            ( "details", "Diagnostics" ),
        )
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GEOS LLVM source coverage</title>
<style>
:root {{ color-scheme: light; --ink: #172033; --muted: #687386; --line: #dfe5ee; --panel: #ffffff; --canvas: #f5f7fb; --blue: #2563eb; --green: #15803d; --red: #b42318; --amber: #b45309; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--canvas); color: var(--ink); font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
main {{ max-width: 1240px; margin: 0 auto; padding: 28px 22px 64px; }}
.hero {{ color: white; border-radius: 16px; padding: 28px 32px; background: linear-gradient(135deg, #172554, #1d4ed8); box-shadow: 0 12px 30px #17255424; }}
.hero.fail {{ background: linear-gradient(135deg, #581c1c, #b42318); }}
.hero.info {{ background: linear-gradient(135deg, #334155, #64748b); }}
.eyebrow {{ margin: 0 0 5px; opacity: .78; text-transform: uppercase; letter-spacing: .12em; font-size: 12px; font-weight: 700; }}
h1, h2, h3 {{ line-height: 1.2; }}
h1 {{ margin: 0; font-size: clamp(28px, 4vw, 42px); }}
h2 {{ margin: 34px 0 10px; font-size: 25px; }}
h3 {{ margin: 26px 0 10px; font-size: 18px; }}
.hero-meta {{ margin-top: 13px; opacity: .9; }}
.hero .badge {{ margin-left: 12px; vertical-align: 2px; }}
nav {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 18px 0 2px; }}
nav a {{ color: var(--blue); background: var(--panel); border: 1px solid var(--line); border-radius: 999px; padding: 5px 12px; text-decoration: none; font-weight: 600; }}
.cards {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin-top: 20px; }}
.card {{ min-height: 144px; background: var(--panel); border: 1px solid var(--line); border-top: 4px solid var(--blue); border-radius: 12px; padding: 17px; box-shadow: 0 3px 10px #1725540b; }}
.card.pass {{ border-top-color: var(--green); }} .card.fail {{ border-top-color: var(--red); }} .card.warn {{ border-top-color: var(--amber); }}
.card-title {{ color: var(--muted); font-weight: 700; font-size: 13px; text-transform: uppercase; letter-spacing: .04em; }}
.card-result {{ margin: 12px 0 8px; font-size: 19px; font-weight: 700; }}
.card-detail, .section-intro, .caption {{ color: var(--muted); }}
.card-detail {{ font-size: 13px; }}
.badge {{ display: inline-flex; align-items: center; gap: 6px; border-radius: 999px; padding: 3px 9px; font-size: 12px; font-weight: 700; white-space: nowrap; }}
.badge-icon {{ display: inline-grid; place-items: center; width: 16px; height: 16px; border-radius: 50%; color: white; background: currentColor; font-size: 11px; }}
.badge.pass {{ color: var(--green); background: #dcfce7; }} .badge.pass .badge-icon {{ color: white; background: var(--green); }}
.badge.fail {{ color: var(--red); background: #fee2e2; }} .badge.fail .badge-icon {{ color: white; background: var(--red); }}
.badge.warn {{ color: var(--amber); background: #fef3c7; }} .badge.warn .badge-icon {{ color: white; background: var(--amber); }}
.badge.info {{ color: #475569; background: #e2e8f0; }} .badge.info .badge-icon {{ color: white; background: #64748b; }}
section {{ background: var(--panel); border: 1px solid var(--line); border-radius: 12px; padding: 22px; margin-top: 18px; box-shadow: 0 3px 10px #1725540b; }}
section h2 {{ margin-top: 0; }}
code {{ color: #334155; background: #eef2f7; border-radius: 4px; padding: 1px 4px; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .9em; }}
.table-wrap {{ overflow-x: auto; }}
table {{ width: 100%; border-collapse: collapse; min-width: 620px; }}
th, td {{ border-bottom: 1px solid var(--line); padding: 10px 9px; text-align: left; vertical-align: middle; }}
th {{ color: var(--muted); background: #f8fafc; font-size: 12px; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }}
tbody tr:last-child td {{ border-bottom: 0; }}
.metrics td:first-child {{ font-weight: 700; }}
.progress {{ width: 110px; height: 6px; margin-top: 5px; overflow: hidden; border-radius: 99px; background: #e2e8f0; }}
.progress span {{ display: block; height: 100%; border-radius: inherit; background: linear-gradient(90deg, #60a5fa, var(--blue)); }}
.delta.pass {{ color: var(--green); font-weight: 700; }} .delta.warn {{ color: var(--amber); font-weight: 700; }}
.caption {{ margin: 12px 0 0; font-size: 13px; }}
.diagnostic-callout {{ margin-top: 14px; padding: 11px 13px; border-left: 4px solid var(--blue); background: #eff6ff; color: #334155; }}
.empty {{ color: var(--muted); padding: 12px 0; }}
.config {{ max-width: 760px; }}
details {{ border-top: 1px solid var(--line); padding: 12px 0; }}
details:first-of-type {{ border-top: 0; }}
summary {{ cursor: pointer; font-weight: 700; }}
.file-name, .file-size {{ color: var(--muted); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; font-weight: 400; margin-left: 8px; }}
pre {{ max-height: 520px; overflow: auto; margin: 12px 0 0; padding: 14px; border-radius: 8px; background: #0f172a; color: #dbeafe; font: 12px/1.45 ui-monospace, SFMono-Regular, Menlo, monospace; white-space: pre-wrap; overflow-wrap: anywhere; }}
@media (max-width: 980px) {{ .cards {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }} }}
@media (max-width: 560px) {{ main {{ padding: 16px 10px 42px; }} .hero {{ padding: 22px; }} .cards {{ grid-template-columns: 1fr; }} section {{ padding: 16px; }} }}
</style>
</head>
<body>
<main>
<header class="hero {escape( overall_state )}">
  <p class="eyebrow">GEOS CI · LLVM source coverage</p>
  <h1>{badge( overall_state, overall_result )} Coverage report</h1>
  <div class="hero-meta">Production scope <code>{scope}</code> · measured commit <code>{commit}</code></div>
</header>
<nav>{nav}</nav>
<div class="cards">
{card( "Coverage smoke tests", test_result, escape( test_detail ), "pass" if "pass" in test_result.lower() else "fail" if "fail" in test_result.lower() else "info" )}
{card( "Enforced policy", policy_result_text, escape( policy_detail ), policy_state )}
{card( "Baseline comparison", baseline_result, escape( baseline_detail ), baseline_state )}
{card( "Changed production lines", changed_result, escape( changed_detail ), changed_state )}
</div>
{coverage_section( summary, thresholds )}
{pull_request_section( report )}
<section id="configuration"><h2>Run configuration</h2>{table( [ "Item", "Value" ], config_rows, "config" ) if config_rows else '<p class="empty">Run configuration is unavailable.</p>'}</section>
{details_section( artifact_dir )}
</main>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser( description=__doc__ )
    parser.add_argument( "--artifact-dir", required=True, type=pathlib.Path )
    parser.add_argument( "--output", required=True, type=pathlib.Path )
    arguments = parser.parse_args()
    try:
        write_atomic( arguments.output, render_html( arguments.artifact_dir ) )
    except ( OSError, ValueError ) as error:
        print( f"Coverage HTML report error: {error}", file=sys.stderr )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit( main() )
