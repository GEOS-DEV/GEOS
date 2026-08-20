#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

import importlib.util
import pathlib
import unittest


SCRIPT = pathlib.Path( __file__ ).parents[1] / "render_coverage_review_summary.py"
SPEC = importlib.util.spec_from_file_location( "coverage_review_summary", SCRIPT )
coverage_review_summary = importlib.util.module_from_spec( SPEC )
SPEC.loader.exec_module( coverage_review_summary )


def metric( covered: int, total: int ) -> dict:
    return {
        "covered": covered,
        "not_covered": total - covered,
        "total": total,
        "percent": round( 100.0 * covered / total, 6 ),
    }


def summary() -> dict:
    return {
        "scope": "src/coreComponents",
        "metrics": {
            "regions": metric( 50, 100 ),
            "functions": metric( 76, 100 ),
            "lines": metric( 74, 100 ),
            "branches": metric( 55, 100 ),
        },
    }


def unavailable_report() -> dict:
    return {
        "aggregate_comparison": {
            "comparable": False,
            "reason": "No exact-base coverage summary was provided.",
        },
        "patch": {
            "metrics": { "executable_lines": metric( 0, 1 ) },
            "top_impacted_locations": [
                {
                    "path": "src/coreComponents/example.cpp",
                    "line": 42,
                    "executable_line": True,
                    "line_covered": False,
                }
            ],
        },
    }


class CoverageReviewSummaryTests( unittest.TestCase ):
    thresholds = { "lines": 7001, "functions": 7501, "branches": 5001 }

    def test_console_summary_is_plain_text_and_contains_enforced_metrics( self ):
        rendered = coverage_review_summary.render_console_summary(
            summary(), self.thresholds
        )

        self.assertIn( "LLVM source coverage policy", rendered )
        self.assertIn( "Overall result: ✅ PASS", rendered )
        self.assertIn( "Lines", rendered )
        self.assertIn( "74.00%", rendered )
        self.assertIn( "74 / 100", rendered )
        self.assertIn( ">= 70.01%", rendered )
        self.assertIn( "Functions", rendered )
        self.assertIn( "Canonical branches", rendered )
        self.assertNotIn( "|---", rendered )

    def test_console_summary_marks_a_failed_threshold( self ):
        failing_summary = summary()
        failing_summary["metrics"]["lines"]["covered"] = 69
        failing_summary["metrics"]["lines"]["not_covered"] = 31
        failing_summary["metrics"]["lines"]["percent"] = 69.0

        rendered = coverage_review_summary.render_console_summary(
            failing_summary, self.thresholds
        )

        self.assertIn( "Overall result: ❌ FAIL", rendered )
        self.assertIn( "❌ Fail", rendered )

    def test_summary_separates_policy_baseline_and_changed_line_signals( self ):
        rendered = coverage_review_summary.render_summary(
            summary(),
            self.thresholds,
            "✅ Pass",
            "106/106 passed in 106.37 s",
            unavailable_report(),
        )

        self.assertIn( "### Reviewer summary", rendered )
        self.assertIn( "106/106 passed in 106.37 s", rendered )
        self.assertIn( "✅ Lines 74.00% (74/100) ≥ 70.01%", rendered )
        self.assertIn( "✅ Functions 76.00% (76/100) ≥ 75.01%", rendered )
        self.assertIn( "✅ Canonical branches 55.00% (55/100) ≥ 50.01%", rendered )
        self.assertIn( "ℹ️ Unavailable", rendered )
        self.assertIn( "No exact-base baseline artifact was available.", rendered )
        self.assertIn( "⚠️ Review", rendered )
        self.assertIn( "example.cpp:42", rendered )

    def test_summary_reports_a_project_coverage_regression( self ):
        report = unavailable_report()
        report["aggregate_comparison"] = {
            "comparable": True,
            "metrics": {
                name: { "percentage_point_delta": delta }
                for name, delta in (
                    ( "lines", -0.10 ),
                    ( "functions", -0.05 ),
                    ( "branches", -0.20 ),
                    ( "regions", -0.01 ),
                )
            },
        }

        rendered = coverage_review_summary.render_summary(
            summary(), self.thresholds, "✅ Pass", "106/106 passed", report
        )

        self.assertIn( "⚠️ Decreased", rendered )
        self.assertIn( "Lines -0.10 pp", rendered )
        self.assertIn( "Canonical branches -0.20 pp", rendered )

    def test_non_pull_request_summary_has_no_pr_warning( self ):
        rendered = coverage_review_summary.render_summary(
            summary(), self.thresholds, "✅ Pass", "106/106 passed", None
        )

        self.assertIn( "⏭️ Not applicable", rendered )
        self.assertIn( "not a pull-request coverage run", rendered )


if __name__ == "__main__":
    unittest.main()
