#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

import importlib.util
import json
import pathlib
import tempfile
import unittest


SCRIPT = pathlib.Path( __file__ ).parents[1] / "render_coverage_report_html.py"
SPEC = importlib.util.spec_from_file_location( "coverage_html", SCRIPT )
coverage_html = importlib.util.module_from_spec( SPEC )
SPEC.loader.exec_module( coverage_html )


def metric( covered: int, total: int ) -> dict:
    return {
        "covered": covered,
        "not_covered": total - covered,
        "total": total,
        "percent": round( 100.0 * covered / total, 6 ) if total else 0.0,
    }


def summary() -> dict:
    return {
        "schema_version": 3,
        "scope": "src/coreComponents",
        "tool": { "name": "llvm-cov", "major": 20 },
        "metrics": {
            "lines": metric( 74, 100 ),
            "functions": metric( 76, 100 ),
            "branches": metric( 55, 100 ),
            "regions": metric( 50, 100 ),
        },
        "inputs": {
            "profiles": 106,
            "coverage_objects": 29,
            "zero_hash_mappings": 3,
        },
        "supplemental": { "native_branch_outcomes": metric( 7, 100 ) },
        "top_branch_gaps": [
            {
                "path": "src/coreComponents/example.cpp",
                **metric( 1, 2 ),
            }
        ],
        "measurement": { "commit_sha": "a" * 40 },
    }


def pull_request_report() -> dict:
    return {
        "base_sha": "b" * 40,
        "head_sha": "a" * 40,
        "scope": "src/coreComponents",
        "aggregate_comparison": {
            "comparable": False,
            "reason": "No exact-base coverage summary was provided.",
        },
        "patch": {
            "changed_files": 1,
            "changed_lines": 2,
            "non_instrumented_changed_lines": 0,
            "metrics": {
                "executable_lines": metric( 0, 1 ),
                "function_start_sites": metric( 0, 0 ),
                "native_branch_outcomes": metric( 0, 0 ),
            },
            "top_impacted_files": [],
            "top_impacted_locations": [
                {
                    "path": "src/coreComponents/example.cpp",
                    "line": 42,
                    "executable_line": True,
                    "line_covered": False,
                    "native_branch_outcomes": metric( 0, 0 ),
                    "function_start_covered": None,
                    "hint": "Add a focused coverage smoke.",
                }
            ],
        },
    }


class CoverageHtmlTests( unittest.TestCase ):
    def test_report_is_a_single_self_contained_human_readable_page( self ):
        with tempfile.TemporaryDirectory() as temporary:
            artifact_dir = pathlib.Path( temporary )
            ( artifact_dir / "coverage-summary.json" ).write_text(
                json.dumps( summary() ), encoding="utf-8"
            )
            ( artifact_dir / "pr-coverage.json" ).write_text(
                json.dumps( pull_request_report() ), encoding="utf-8"
            )
            ( artifact_dir / "coverage-status.md" ).write_text(
                """### ✅ LLVM source coverage passed

| Validation | Result | Details |
|---|:---:|---|
| Coverage smoke suite | ✅ Pass | 106/106 passed in 10.00 s |

### Run configuration

| Item | Value |
|---|---|
| Toolchain | `Clang 20 / LLVM 20` |
""",
                encoding="utf-8",
            )
            ( artifact_dir / "ctest.log" ).write_text(
                "100% tests passed, 0 tests failed out of 106\n", encoding="utf-8"
            )
            ( artifact_dir / "mapping-integrity.log" ).write_text(
                "mapping checks passed\n", encoding="utf-8"
            )
            output = coverage_html.render_html( artifact_dir )

        self.assertIn( "<!doctype html>", output )
        self.assertIn( "Enforced policy", output )
        self.assertIn( "74.00%", output )
        self.assertIn( "No exact-base coverage summary was provided.", output )
        self.assertIn( "example.cpp:42", output )
        self.assertIn( "mapping checks passed", output )
        self.assertIn( "details", output )


if __name__ == "__main__":
    unittest.main()
