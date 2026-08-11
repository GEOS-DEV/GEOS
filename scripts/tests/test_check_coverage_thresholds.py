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
        "percent": round( 100.0 * covered / total if total else 0.0, 6 ),
    }


def measurement( scope, excluded_regex ):
    value = {
        "commit_sha": "1" * 40,
        "tree_sha": "2" * 40,
        "contract_id": coverage_gate.COVERAGE_CONTRACT_ID,
        "contract_fingerprint": "0" * 64,
        "container": {
            "image": "geosx/ubuntu24.04-clang20:349-1014",
            "image_id": f"sha256:{'3' * 64}",
            "image_digests": [
                f"geosx/ubuntu24.04-clang20@sha256:{'4' * 64}"
            ],
        },
        "toolchain": {
            "c_compiler_version": (
                "Ubuntu clang version 20.1.2 (0ubuntu1~24.04.3)\n"
                "Target: x86_64-pc-linux-gnu"
            ),
            "cxx_compiler_version": (
                "Ubuntu clang version 20.1.2 (0ubuntu1~24.04.3)\n"
                "Target: x86_64-pc-linux-gnu"
            ),
            "llvm_cov_version": (
                "Ubuntu LLVM version 20.1.2\n  Optimized build."
            ),
            "llvm_package_versions": [
                "libclang-rt-20-dev=1:20.1.2-0ubuntu1~24.04.3",
                "llvm-20=1:20.1.2-0ubuntu1~24.04.3",
            ],
            "compiler_target": "x86_64-pc-linux-gnu",
        },
        "build_config": {
            "BLT_CXX_STD": "c++17",
            "BUILD_SHARED_LIBS": True,
            "CMAKE_BUILD_TYPE": "RelWithDebInfo",
            "CMAKE_C_FLAGS": "-pthread",
            "CMAKE_C_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
            "CMAKE_CXX_FLAGS": "-pthread",
            "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
            "ENABLE_COVERAGE": False,
            "ENABLE_CUDA": False,
            "ENABLE_HIP": False,
            "ENABLE_HYPRE": True,
            "ENABLE_MPI": True,
            "ENABLE_OPENMP": True,
            "ENABLE_HYPRE_DEVICE": "CPU",
            "ENABLE_TRILINOS": False,
            "GEOS_BUILD_SHARED_LIBS": True,
            "GEOS_ENABLE_BOUNDS_CHECK": False,
            "GEOS_ENABLE_LLVM_SOURCE_COVERAGE": True,
            "GEOS_GLOBALINDEX_TYPE": "long long int",
            "GEOS_LA_INTERFACE": "Hypre",
            "GEOS_LOCALINDEX_TYPE": "int",
            "LVARRAY_BOUNDS_CHECK": False,
            "RAJA_ENABLE_CUDA": False,
            "RAJA_ENABLE_HIP": False,
            "RAJA_ENABLE_OPENMP": False,
        },
        "metric_semantics": json.loads(
            json.dumps( coverage_gate.METRIC_SEMANTICS )
        ),
        "object_selection": json.loads(
            json.dumps( coverage_gate.OBJECT_SELECTION )
        ),
    }
    value["contract_fingerprint"] = coverage_gate.compute_contract_fingerprint(
        scope, excluded_regex, value
    )
    return value


def summary( branch_covered=1, branch_total=2 ):
    scope = "src/coreComponents"
    excluded_regex = "/unitTests/"
    canonical = {
        "regions": metric( 1, 2 ),
        "functions": metric( 1, 2 ),
        "lines": metric( 1, 2 ),
        "branches": metric( branch_covered, branch_total ),
    }
    document = {
        "schema_version": 3,
        "scope": scope,
        "excluded_regex": excluded_regex,
        "measurement": measurement( scope, excluded_regex ),
        "tool": { "name": "llvm-cov", "major": 20 },
        "inputs": {
            "profiles": 10,
            "coverage_objects": 2,
            "zero_hash_mappings": 3,
        },
        "metrics": canonical,
        "per_file_metrics": [
            {
                "path": "src/coreComponents/example.cpp",
                "metrics": json.loads( json.dumps( canonical ) ),
            }
        ],
        "top_branch_gaps": [],
        "supplemental": {
            "native_branch_outcomes": metric( 1, 4 ),
        },
    }
    if branch_covered < branch_total:
        document["top_branch_gaps"] = [
            {
                "path": "src/coreComponents/example.cpp",
                **metric( branch_covered, branch_total ),
            }
        ]
    return document


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
        with self.assertRaisesRegex( ValueError, "integer 3" ):
            coverage_gate.validate_summary( document )

    def test_invalid_metadata_is_rejected( self ):
        document = summary()
        document["inputs"]["profiles"] = 0
        with self.assertRaisesRegex( ValueError, "inputs.profiles" ):
            coverage_gate.validate_summary( document )

    def test_contract_fingerprint_is_stable_and_excludes_revision_inputs( self ):
        document = summary()
        self.assertEqual(
            document["measurement"]["contract_fingerprint"],
            "392f58cf380615f8214dc716273efbadcabefb0f25d6bfdd83560077f1cd4672",
        )

        document["measurement"]["commit_sha"] = "a" * 40
        document["measurement"]["tree_sha"] = "b" * 40
        document["measurement"]["container"]["image"] = "local/relabel:latest"
        document["inputs"]["profiles"] = 999
        document["inputs"]["coverage_objects"] = 77
        validated = coverage_gate.validate_summary( document )
        self.assertEqual( validated["measurement"]["commit_sha"], "a" * 40 )

    def test_hard_contract_changes_require_a_new_fingerprint( self ):
        mutations = (
            lambda document: document["measurement"]["build_config"].__setitem__(
                "ENABLE_HYPRE", False
            ),
            lambda document: document["measurement"]["toolchain"].__setitem__(
                "compiler_target", "aarch64-unknown-linux-gnu"
            ),
            lambda document: document["measurement"]["toolchain"].__setitem__(
                "llvm_package_versions", [ "llvm-20=20.1.3" ]
            ),
            lambda document: document["measurement"]["container"].__setitem__(
                "image_id", f"sha256:{'5' * 64}"
            ),
        )
        for mutate in mutations:
            with self.subTest( mutation=mutate ):
                document = summary()
                mutate( document )
                with self.assertRaisesRegex( ValueError, "fingerprint" ):
                    coverage_gate.validate_summary( document )

    def test_measurement_contract_fields_are_strictly_validated( self ):
        document = summary()
        document["measurement"]["contract_id"] = "other-contract"
        with self.assertRaisesRegex( ValueError, "contract_id" ):
            coverage_gate.validate_summary( document )

        document = summary()
        document["measurement"]["toolchain"]["llvm_cov_version"] = (
            "Ubuntu LLVM version 19.1.0"
        )
        with self.assertRaisesRegex( ValueError, "does not match tool.major" ):
            coverage_gate.validate_summary( document )

        document = summary()
        package = document["measurement"]["toolchain"]["llvm_package_versions"][0]
        document["measurement"]["toolchain"]["llvm_package_versions"] = [
            package,
            package,
        ]
        with self.assertRaisesRegex( ValueError, "sorted and unique" ):
            coverage_gate.validate_summary( document )

        document = summary()
        digest = document["measurement"]["container"]["image_digests"][0]
        document["measurement"]["container"]["image_digests"] = [ digest, digest ]
        with self.assertRaisesRegex( ValueError, "sorted and unique" ):
            coverage_gate.validate_summary( document )

    def test_build_config_must_be_normalized_and_relevant( self ):
        document = summary()
        document["measurement"]["build_config"]["ENABLE_COVERAGE"] = "OFF"
        with self.assertRaisesRegex( ValueError, "ENABLE_COVERAGE" ):
            coverage_gate.validate_summary( document )

        document = summary()
        document["measurement"]["build_config"]["CMAKE_INSTALL_PREFIX"] = "/tmp"
        with self.assertRaisesRegex( ValueError, "irrelevant key" ):
            coverage_gate.validate_summary( document )

    def test_inconsistent_percent_is_rejected( self ):
        document = summary()
        document["metrics"]["branches"]["percent"] = 99.0
        with self.assertRaisesRegex( ValueError, "percent is inconsistent" ):
            coverage_gate.validate_summary( document )

    def test_per_file_metrics_are_sorted_unique_and_in_scope( self ):
        document = summary()
        document["per_file_metrics"] = [
            {
                "path": "src/coreComponents/z.cpp",
                "metrics": {
                    name: metric( 0, 1 ) for name in coverage_gate.CANONICAL_METRICS
                },
            },
            {
                "path": "src/coreComponents/a.cpp",
                "metrics": {
                    name: metric( 1, 1 ) for name in coverage_gate.CANONICAL_METRICS
                },
            },
        ]
        document["top_branch_gaps"] = [
            { "path": "src/coreComponents/z.cpp", **metric( 0, 1 ) }
        ]
        with self.assertRaisesRegex( ValueError, "sorted by path" ):
            coverage_gate.validate_summary( document )

        document = summary()
        document["per_file_metrics"][0]["path"] = "outside.cpp"
        with self.assertRaisesRegex( ValueError, "outside the scope" ):
            coverage_gate.validate_summary( document )

        document = summary()
        duplicate = json.loads( json.dumps( document["per_file_metrics"][0] ) )
        document["per_file_metrics"].append( duplicate )
        with self.assertRaisesRegex( ValueError, "paths must be unique" ):
            coverage_gate.validate_summary( document )

    def test_per_file_metrics_must_reconcile_every_canonical_metric( self ):
        for name in coverage_gate.CANONICAL_METRICS:
            with self.subTest( metric=name ):
                document = summary()
                document["per_file_metrics"][0]["metrics"][name] = metric( 0, 2 )
                with self.assertRaisesRegex( ValueError, f"{name} covered" ):
                    coverage_gate.validate_summary( document )

    def test_per_file_zero_denominators_are_allowed( self ):
        document = summary()
        document["per_file_metrics"].insert(
            0,
            {
                "path": "src/coreComponents/empty.cpp",
                "metrics": {
                    name: metric( 0, 0 ) for name in coverage_gate.CANONICAL_METRICS
                },
            },
        )
        validated = coverage_gate.validate_summary( document )
        self.assertEqual(
            validated["per_file_metrics"][0]["metrics"]["branches"]["total"],
            0,
        )

    def test_top_branch_gaps_must_match_per_file_metrics( self ):
        document = summary()
        document["top_branch_gaps"][0] = {
            "path": "src/coreComponents/example.cpp",
            **metric( 0, 2 ),
        }
        with self.assertRaisesRegex( ValueError, "does not match canonical per-file" ):
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
| <code>example.cpp</code> | 1 | 1 | 2 | 50.00% |

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
        document = summary( 12, 23 )
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
                **metric( 10, 20 ),
            },
            {
                "path": "src/coreComponents/nested/b.cpp",
                **metric( 2, 3 ),
            },
        ]
        document["per_file_metrics"] = [
            {
                "path": "src/coreComponents/a.cpp",
                "metrics": {
                    "regions": metric( 1, 2 ),
                    "functions": metric( 10887, 14276 ),
                    "lines": metric( 1, 2 ),
                    "branches": metric( 10, 20 ),
                },
            },
            {
                "path": "src/coreComponents/nested/b.cpp",
                "metrics": {
                    "regions": metric( 0, 0 ),
                    "functions": metric( 0, 0 ),
                    "lines": metric( 0, 0 ),
                    "branches": metric( 2, 3 ),
                },
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
        self.assertIn( "| <code>a.cpp</code> | 10 | 10 | 20 | 50.00% |", markdown )
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
                with self.assertRaisesRegex( ValueError, "normalized" ):
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
