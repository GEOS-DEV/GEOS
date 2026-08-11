#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

import importlib.util
import json
import pathlib
import unittest


SCRIPT = pathlib.Path( __file__ ).parents[1] / "check_coverage_thresholds.py"
REPOSITORY_ROOT = pathlib.Path( __file__ ).parents[2]
SPEC = importlib.util.spec_from_file_location( "coverage_gate", SCRIPT )
coverage_gate = importlib.util.module_from_spec( SPEC )
SPEC.loader.exec_module( coverage_gate )


def metric( covered, total ):
    return {
        "covered": covered,
        "total": total,
        "not_covered": total - covered,
        "percent": round( 100.0 * covered / total, 6 ),
    }


def summary( branch_covered=1, branch_total=2 ):
    return {
        "schema_version": 2,
        "scope": "src/coreComponents",
        "excluded_regex": "/unitTests/",
        "tool": { "name": "llvm-cov", "major": 20 },
        "inputs": {
            "profiles": 10,
            "coverage_objects": 2,
            "zero_hash_mappings": 3,
        },
        "metrics": {
            "regions": metric( 1, 2 ),
            "functions": metric( 1, 2 ),
            "lines": metric( 1, 2 ),
            "branches": metric( branch_covered, branch_total ),
        },
        "top_branch_gaps": [
            {
                "path": "src/coreComponents/example.cpp",
                **metric( 0, 1 ),
            }
        ],
        "supplemental": {
            "native_branch_outcomes": metric( 1, 4 ),
        },
    }


class CoverageThresholdTests( unittest.TestCase ):
    def test_exact_threshold_passes( self ):
        validated = coverage_gate.validate_summary( summary() )
        markdown, passed = coverage_gate.build_markdown(
            validated, { "branches": 5000 }
        )
        self.assertTrue( passed )
        self.assertIn( "50.00%", markdown )
        self.assertIn( "✅ Coverage policy passed", markdown )

    def test_one_count_below_threshold_fails( self ):
        validated = coverage_gate.validate_summary( summary( 49, 100 ) )
        markdown, passed = coverage_gate.build_markdown(
            validated, { "branches": 5000 }
        )
        self.assertFalse( passed )
        self.assertIn( "❌ Coverage policy failed", markdown )
        self.assertIn( "❌ Fail (-1.00 pp)", markdown )

    def test_integer_comparison_avoids_rounding_up( self ):
        validated = coverage_gate.validate_summary( summary( 1, 3 ) )
        markdown, passed = coverage_gate.build_markdown(
            validated, { "branches": 3334 }
        )
        self.assertFalse( passed )
        self.assertIn( "❌ Fail (<0.01 pp below)", markdown )

    def test_checked_in_threshold_boundary_is_exact( self ):
        passing = coverage_gate.validate_summary( summary( 5001, 10000 ) )
        failing = coverage_gate.validate_summary( summary( 5000, 10000 ) )
        self.assertTrue(
            coverage_gate.build_markdown( passing, { "branches": 5001 } )[1]
        )
        self.assertFalse(
            coverage_gate.build_markdown( failing, { "branches": 5001 } )[1]
        )

    def test_inconsistent_counts_are_rejected( self ):
        document = summary()
        document["metrics"]["branches"]["not_covered"] = 99
        with self.assertRaisesRegex( ValueError, "inconsistent" ):
            coverage_gate.validate_summary( document )

    def test_boolean_schema_version_is_rejected( self ):
        document = summary()
        document["schema_version"] = True
        with self.assertRaisesRegex( ValueError, "integer 2" ):
            coverage_gate.validate_summary( document )

    def test_invalid_metadata_is_rejected( self ):
        document = summary()
        document["inputs"]["profiles"] = 0
        with self.assertRaisesRegex( ValueError, "inputs.profiles" ):
            coverage_gate.validate_summary( document )

    def test_inconsistent_percent_is_rejected( self ):
        document = summary()
        document["metrics"]["branches"]["percent"] = 99.0
        with self.assertRaisesRegex( ValueError, "percent is inconsistent" ):
            coverage_gate.validate_summary( document )

    def test_policy_scope_must_match( self ):
        policy = {
            "schema_version": 1,
            "scope": "wrong/scope",
            "minimum_basis_points": { "branches": 5000 },
        }
        with self.assertRaisesRegex( ValueError, "does not match" ):
            coverage_gate.validate_policy( policy, "src/coreComponents" )

    def test_non_gated_metrics_are_report_only( self ):
        validated = coverage_gate.validate_summary( summary() )
        markdown, passed = coverage_gate.build_markdown(
            validated, { "branches": 5000 }
        )
        self.assertTrue( passed )
        self.assertIn( "Functions", markdown )
        self.assertIn( "ℹ️ Measured", markdown )

    def test_markdown_golden_output( self ):
        validated = coverage_gate.validate_summary( summary() )
        markdown, passed = coverage_gate.build_markdown(
            validated, { "branches": 5000 }
        )
        expected = """### ✅ Coverage policy passed

Production scope: `src/coreComponents` · Reporter: `llvm-cov 20`

| Metric | Covered | Missed | Total | Coverage | Requirement | Result |
|---|---:|---:|---:|---:|---:|:---:|
| Lines | 1 | 1 | 2 | 50.00% | — | ℹ️ Measured |
| Functions | 1 | 1 | 2 | 50.00% | — | ℹ️ Measured |
| Canonical branches | 1 | 1 | 2 | 50.00% | ≥ 50.00% | ✅ Pass (+0.00 pp) |
| Regions | 1 | 1 | 2 | 50.00% | — | ℹ️ Measured |

Percentages are rounded for display; policy thresholds use exact covered/total counts.

### Largest branch gaps

| Source file | Covered | Missed | Total | Coverage |
|---|---:|---:|---:|---:|
| <code>example.cpp</code> | 0 | 1 | 1 | 0.00% |

Ranked by missed canonical branch outcomes within the production scope.

### Coverage collection

| Input | Value |
|---|---:|
| Profile files merged | 10 |
| Instrumented product objects | 2 |
| Audited zero-hash mappings | 3 |

Zero-hash mappings have no profile counters and are audited during report generation; any nonzero or unaudited mismatch fails the job.

<details>
<summary>Branch diagnostic (not gated)</summary>

| Diagnostic | Covered | Total | Coverage |
|---|---:|---:|---:|
| Native branch outcomes | 1 | 4 | 25.00% |

Native outcomes are instantiation-weighted: C++ template and inline-function instances may appear more than once. This diagnostic is not the policy metric.

</details>
"""
        self.assertTrue( passed )
        self.assertEqual( markdown, expected )

    def test_markdown_presents_counts_inputs_and_diagnostics( self ):
        document = summary( 19686, 35810 )
        document["inputs"] = {
            "profiles": 3451,
            "coverage_objects": 29,
            "zero_hash_mappings": 29966,
        }
        document["metrics"]["functions"] = metric( 10887, 14276 )
        document["supplemental"]["native_branch_outcomes"] = metric(
            48401, 684500
        )
        document["top_branch_gaps"] = [
            {
                "path": "src/coreComponents/a.cpp",
                **metric( 10, 30 ),
            },
            {
                "path": "src/coreComponents/nested/b.cpp",
                **metric( 2, 7 ),
            },
        ]
        validated = coverage_gate.validate_summary( document )
        markdown, passed = coverage_gate.build_markdown(
            validated, { "branches": 5001 }
        )

        self.assertTrue( passed )
        self.assertIn(
            "| Functions | 10,887 | 3,389 | 14,276 | 76.26% | — | ℹ️ Measured |",
            markdown,
        )
        self.assertIn( "### Largest branch gaps", markdown )
        self.assertLess( markdown.index( "a.cpp" ), markdown.index( "nested/b.cpp" ) )
        self.assertIn( "| <code>a.cpp</code> | 10 | 20 | 30 | 33.33% |", markdown )
        self.assertIn( "| Profile files merged | 3,451 |", markdown )
        self.assertIn( "| Instrumented product objects | 29 |", markdown )
        self.assertIn( "| Audited zero-hash mappings | 29,966 |", markdown )
        self.assertIn( "| Native branch outcomes | 48,401 | 684,500 | 7.07% |", markdown )
        primary_table, diagnostic = markdown.split( "<details>", maxsplit=1 )
        self.assertNotIn( "Native branch outcomes", primary_table )
        self.assertIn( "instantiation-weighted", diagnostic )

    def test_branch_gap_paths_must_be_unique_and_in_scope( self ):
        document = summary()
        document["top_branch_gaps"][0]["path"] = "src/tests/not-production.cpp"
        with self.assertRaisesRegex( ValueError, "outside the scope" ):
            coverage_gate.validate_summary( document )

        document = summary()
        document["top_branch_gaps"].append( dict( document["top_branch_gaps"][0] ) )
        with self.assertRaisesRegex( ValueError, "normalized and unique" ):
            coverage_gate.validate_summary( document )

        for invalid_path in (
            "src/coreComponents/",
            "src/coreComponents//file.cpp",
            "src/coreComponents/./file.cpp",
        ):
            with self.subTest( path=invalid_path ):
                document = summary()
                document["top_branch_gaps"][0]["path"] = invalid_path
                with self.assertRaisesRegex( ValueError, "normalized and unique" ):
                    coverage_gate.validate_summary( document )

    def test_branch_gaps_must_be_sorted_by_missed_count( self ):
        document = summary()
        document["top_branch_gaps"] = [
            { "path": "src/coreComponents/small.cpp", **metric( 2, 3 ) },
            { "path": "src/coreComponents/large.cpp", **metric( 1, 4 ) },
        ]
        with self.assertRaisesRegex( ValueError, "sorted" ):
            coverage_gate.validate_summary( document )

    def test_branch_gap_totals_cannot_exceed_canonical_branches( self ):
        document = summary()
        document["top_branch_gaps"] = [
            { "path": "src/coreComponents/too-large.cpp", **metric( 0, 3 ) }
        ]
        with self.assertRaisesRegex( ValueError, "exceed canonical" ):
            coverage_gate.validate_summary( document )

    def test_checked_in_policy_is_valid( self ):
        policy_path = REPOSITORY_ROOT / ".github" / "coverage-thresholds.json"
        with policy_path.open( encoding="utf-8" ) as stream:
            policy = json.load( stream )
        thresholds = coverage_gate.validate_policy( policy, "src/coreComponents" )
        self.assertEqual( thresholds, { "branches": 5001 } )


if __name__ == "__main__":
    unittest.main()
