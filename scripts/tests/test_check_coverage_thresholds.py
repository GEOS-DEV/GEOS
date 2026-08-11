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
        "schema_version": 1,
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
        self.assertIn( "50.000000%", markdown )

    def test_one_count_below_threshold_fails( self ):
        validated = coverage_gate.validate_summary( summary( 49, 100 ) )
        _, passed = coverage_gate.build_markdown(
            validated, { "branches": 5000 }
        )
        self.assertFalse( passed )

    def test_integer_comparison_avoids_rounding_up( self ):
        validated = coverage_gate.validate_summary( summary( 1, 3 ) )
        _, passed = coverage_gate.build_markdown(
            validated, { "branches": 3334 }
        )
        self.assertFalse( passed )

    def test_inconsistent_counts_are_rejected( self ):
        document = summary()
        document["metrics"]["branches"]["not_covered"] = 99
        with self.assertRaisesRegex( ValueError, "inconsistent" ):
            coverage_gate.validate_summary( document )

    def test_boolean_schema_version_is_rejected( self ):
        document = summary()
        document["schema_version"] = True
        with self.assertRaisesRegex( ValueError, "integer 1" ):
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
        self.assertIn( "Report only", markdown )

    def test_checked_in_policy_is_valid( self ):
        policy_path = REPOSITORY_ROOT / ".github" / "coverage-thresholds.json"
        with policy_path.open( encoding="utf-8" ) as stream:
            policy = json.load( stream )
        thresholds = coverage_gate.validate_policy( policy, "src/coreComponents" )
        self.assertEqual( thresholds, { "branches": 5001 } )


if __name__ == "__main__":
    unittest.main()
