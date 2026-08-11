#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

from __future__ import annotations

import importlib.util
import json
import pathlib
import stat
import subprocess
import sys
import tempfile
import unittest


SCRIPT = pathlib.Path( __file__ ).parents[1] / "compare_pr_coverage.py"
SPEC = importlib.util.spec_from_file_location( "compare_pr_coverage", SCRIPT )
compare_pr_coverage = importlib.util.module_from_spec( SPEC )
SPEC.loader.exec_module( compare_pr_coverage )


def run_command( repository: pathlib.Path, *arguments: str ) -> str:
    result = subprocess.run(
        [ *arguments ],
        cwd=repository,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def metric( covered: int, total: int ) -> dict:
    return {
        "covered": covered,
        "not_covered": total - covered,
        "total": total,
        "percent": round( 100.0 * covered / total, 6 ),
    }


class GitFixture:
    def __init__( self ) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.repository = pathlib.Path( self.temporary.name ) / "repository"
        self.repository.mkdir()
        run_command( self.repository, "git", "init", "-q", "-b", "develop" )
        run_command(
            self.repository, "git", "config", "user.email", "coverage@example.com"
        )
        run_command(
            self.repository, "git", "config", "user.name", "Coverage Tests"
        )

    def close( self ) -> None:
        self.temporary.cleanup()

    def write( self, path: str, contents: str ) -> pathlib.Path:
        target = self.repository / path
        target.parent.mkdir( parents=True, exist_ok=True )
        target.write_text( contents, encoding="utf-8" )
        return target

    def commit( self, message: str ) -> str:
        run_command( self.repository, "git", "add", "--all" )
        run_command( self.repository, "git", "commit", "-q", "-m", message )
        return run_command( self.repository, "git", "rev-parse", "HEAD" )

    def tree( self, commit_sha: str ) -> str:
        return run_command(
            self.repository, "git", "rev-parse", f"{commit_sha}^{{tree}}"
        )


def summary_document(
    commit_sha: str,
    tree_sha: str,
    *,
    schema_version: int = 3,
    image_digit: str = "3",
    scope: str = "src/coreComponents",
    covered_offset: int = 0,
) -> dict:
    document = {
        "schema_version": schema_version,
        "scope": scope,
        "excluded_regex": "/(unitTests|tests)/",
        "tool": { "name": "llvm-cov", "major": 20 },
        "inputs": {
            "profiles": 4,
            "coverage_objects": 2,
            "zero_hash_mappings": 0,
        },
        "metrics": {
            "regions": metric( 50 + covered_offset, 100 ),
            "functions": metric( 60 + covered_offset, 100 ),
            "lines": metric( 70 + covered_offset, 100 ),
            "branches": metric( 40 + covered_offset, 100 ),
        },
        "top_branch_gaps": [],
        "supplemental": {
            "native_branch_outcomes": metric( 5 + covered_offset, 100 )
        },
    }
    if schema_version >= 3:
        measurement = {
            "commit_sha": commit_sha,
            "tree_sha": tree_sha,
            "contract_id": compare_pr_coverage.coverage_contract.COVERAGE_CONTRACT_ID,
            "contract_fingerprint": "0" * 64,
            "container": {
                "image": "geosx/ubuntu24.04-clang20:test",
                "image_id": f"sha256:{image_digit * 64}",
                "image_digests": [
                    f"geosx/ubuntu24.04-clang20@sha256:{'4' * 64}"
                ],
            },
            "toolchain": {
                "c_compiler_version": "Ubuntu clang version 20.1.2",
                "cxx_compiler_version": "Ubuntu clang version 20.1.2",
                "llvm_cov_version": "Ubuntu LLVM version 20.1.2",
                "llvm_package_versions": [
                    "libclang-rt-20-dev=1:20.1.2",
                    "llvm-20=1:20.1.2",
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
                json.dumps( compare_pr_coverage.coverage_contract.METRIC_SEMANTICS )
            ),
            "object_selection": json.loads(
                json.dumps( compare_pr_coverage.coverage_contract.OBJECT_SELECTION )
            ),
        }
        measurement["contract_fingerprint"] = (
            compare_pr_coverage.coverage_contract.compute_contract_fingerprint(
                scope, document["excluded_regex"], measurement
            )
        )
        document["measurement"] = measurement
        document["per_file_metrics"] = [
            {
                "path": "src/coreComponents/New.cpp",
                "metrics": json.loads( json.dumps( document["metrics"] ) ),
            }
        ]
        document["top_branch_gaps"] = [
            {
                "path": "src/coreComponents/New.cpp",
                **document["metrics"]["branches"],
            }
        ]
    return document


class ComparePrCoverageTest( unittest.TestCase ):
    def setUp( self ) -> None:
        self.fixture = GitFixture()

    def tearDown( self ) -> None:
        self.fixture.close()

    def write_json( self, name: str, document: dict ) -> pathlib.Path:
        path = pathlib.Path( self.fixture.temporary.name ) / name
        path.write_text( json.dumps( document ), encoding="utf-8" )
        return path

    def write_lcov( self, contents: str ) -> pathlib.Path:
        path = pathlib.Path( self.fixture.temporary.name ) / "native.info"
        path.write_text( contents, encoding="utf-8" )
        return path

    def make_added_source( self, path: str = "src/coreComponents/New.cpp" ) -> tuple:
        self.fixture.write( "README", "base\n" )
        base_sha = self.fixture.commit( "base" )
        source = self.fixture.write(
            path,
            "int uncovered() {\n"
            "  if ( flag ) {\n"
            "    return 1;\n"
            "  }\n"
            "  return 0;\n"
            "}\n",
        )
        head_sha = self.fixture.commit( "add source" )
        return base_sha, head_sha, source

    def candidate_summary(
        self, head_sha: str, *, schema_version: int = 3, offset: int = 0
    ) -> pathlib.Path:
        return self.write_json(
            "candidate.json",
            summary_document(
                head_sha,
                self.fixture.tree( head_sha ),
                schema_version=schema_version,
                covered_offset=offset,
            ),
        )

    def test_patch_metrics_and_actionable_locations( self ) -> None:
        base_sha, head_sha, source = self.make_added_source()
        lcov = self.write_lcov(
            f"SF:{source}\n"
            "FN:1,uncovered<int>\n"
            "FNDA:0,uncovered<int>\n"
            "DA:1,0\n"
            "DA:2,1\n"
            "DA:3,0\n"
            "DA:5,1\n"
            "BRDA:2,0,0,1\n"
            "BRDA:2,0,1,0\n"
            "end_of_record\n"
        )
        report, markdown = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha ),
        )

        patch = report["patch"]
        self.assertEqual( patch["changed_files"], 1 )
        self.assertEqual( patch["changed_lines"], 6 )
        self.assertEqual(
            patch["metrics"]["executable_lines"], metric( 2, 4 )
        )
        self.assertEqual(
            patch["metrics"]["function_start_sites"], metric( 0, 1 )
        )
        self.assertEqual(
            patch["metrics"]["native_branch_outcomes"], metric( 1, 2 )
        )
        self.assertEqual( patch["non_instrumented_changed_lines"], 2 )
        self.assertIn( "Add cases that exercise the remaining", markdown )
        self.assertIn( "Advisory and instantiation-weighted", markdown )
        self.assertFalse( report["aggregate_comparison"]["comparable"] )
        self.assertIn(
            "No exact-base", report["aggregate_comparison"]["reason"]
        )

    def test_exact_baseline_comparison( self ) -> None:
        base_sha, head_sha, source = self.make_added_source()
        lcov = self.write_lcov(
            f"SF:{source}\nDA:1,1\nend_of_record\n"
        )
        baseline = self.write_json(
            "baseline.json",
            summary_document( base_sha, self.fixture.tree( base_sha ) ),
        )
        candidate = self.candidate_summary( head_sha, offset=5 )
        report, markdown = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            candidate,
            baseline,
        )

        comparison = report["aggregate_comparison"]
        self.assertTrue( comparison["comparable"] )
        self.assertEqual( comparison["metrics"]["branches"]["covered_delta"], 5 )
        self.assertEqual(
            comparison["metrics"]["branches"]["not_covered_delta"], -5
        )
        self.assertEqual(
            comparison["metrics"]["branches"]["percentage_point_delta"], 5.0
        )
        self.assertIn( "| Canonical branches |", markdown )
        self.assertIn( "+5.00 pp", markdown )
        self.assertIn( "##### Canonical movement in changed files", markdown )
        self.assertIn( "| <code>src/coreComponents/New.cpp</code>", markdown )

    def test_incompatible_baseline_is_explicitly_not_comparable( self ) -> None:
        base_sha, head_sha, source = self.make_added_source()
        lcov = self.write_lcov(
            f"SF:{source}\nDA:1,1\nend_of_record\n"
        )
        baseline = self.write_json(
            "baseline.json",
            summary_document(
                base_sha,
                self.fixture.tree( base_sha ),
                image_digit="5",
            ),
        )
        report, markdown = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha ),
            baseline,
        )

        comparison = report["aggregate_comparison"]
        self.assertFalse( comparison["comparable"] )
        self.assertIn( "fingerprints differ", comparison["reason"] )
        self.assertIn( "ℹ️ N/A", markdown )

    def test_tampered_v3_baseline_is_rejected_by_contract_validator( self ) -> None:
        base_sha, head_sha, source = self.make_added_source()
        lcov = self.write_lcov( f"SF:{source}\nDA:1,1\nend_of_record\n" )
        baseline_document = summary_document(
            base_sha, self.fixture.tree( base_sha )
        )
        baseline_document["measurement"]["container"]["image_id"] = (
            f"sha256:{'9' * 64}"
        )
        baseline = self.write_json( "baseline.json", baseline_document )
        report, _ = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha ),
            baseline,
        )

        comparison = report["aggregate_comparison"]
        self.assertFalse( comparison["comparable"] )
        self.assertIn( "baseline summary is invalid", comparison["reason"] )
        self.assertIn( "fingerprint", comparison["reason"] )

    def test_legacy_v2_summaries_make_exact_delta_unavailable( self ) -> None:
        base_sha, head_sha, source = self.make_added_source()
        lcov = self.write_lcov(
            f"SF:{source}\nDA:1,1\nend_of_record\n"
        )
        baseline = self.write_json(
            "baseline.json",
            summary_document(
                base_sha, self.fixture.tree( base_sha ), schema_version=2
            ),
        )
        report, _ = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha, schema_version=2 ),
            baseline,
        )

        comparison = report["aggregate_comparison"]
        self.assertFalse( comparison["comparable"] )
        self.assertIn( "schema-3 measurement provenance", comparison["reason"] )

    def test_rename_uses_only_modified_head_lines( self ) -> None:
        original = "\n".join( f"int value_{index} = {index};" for index in range( 10 ) )
        self.fixture.write( "src/coreComponents/Old.cpp", original + "\n" )
        base_sha = self.fixture.commit( "base" )
        run_command(
            self.fixture.repository,
            "git",
            "mv",
            "src/coreComponents/Old.cpp",
            "src/coreComponents/New.cpp",
        )
        updated = original.replace( "int value_4 = 4;", "int value_4 = 40;" )
        self.fixture.write( "src/coreComponents/New.cpp", updated + "\n" )
        head_sha = self.fixture.commit( "rename and edit" )

        records = compare_pr_coverage.changed_line_ranges(
            self.fixture.repository, base_sha, head_sha
        )
        self.assertEqual( len( records ), 1 )
        self.assertEqual( records[0]["status"], "R" )
        self.assertEqual( records[0]["old_path"], "src/coreComponents/Old.cpp" )
        self.assertEqual( records[0]["path"], "src/coreComponents/New.cpp" )
        self.assertEqual( records[0]["lines"], frozenset( ( 5, ) ) )

    def test_noninstrumented_only_patch_has_empty_executable_metric( self ) -> None:
        source = self.fixture.write(
            "src/coreComponents/Comment.cpp",
            "// old comment\nint stable() { return 1; }\n",
        )
        base_sha = self.fixture.commit( "base" )
        self.fixture.write(
            "src/coreComponents/Comment.cpp",
            "// new comment\nint stable() { return 1; }\n",
        )
        head_sha = self.fixture.commit( "comment only" )
        lcov = self.write_lcov(
            f"SF:{source}\n"
            "FN:2,stable\nFNDA:1,stable\nDA:2,1\nend_of_record\n"
        )
        report, markdown = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha ),
        )

        lines = report["patch"]["metrics"]["executable_lines"]
        self.assertEqual( lines["total"], 0 )
        self.assertIsNone( lines["percent"] )
        self.assertEqual( report["patch"]["non_instrumented_changed_lines"], 1 )
        self.assertIn( "| Executable changed lines | 0 | 0 | 0 | N/A |", markdown )

    def test_test_sources_are_excluded_from_production_patch_totals( self ) -> None:
        self.fixture.write( "README", "base\n" )
        base_sha = self.fixture.commit( "base" )
        self.fixture.write(
            "src/coreComponents/unitTests/testFeature.cpp",
            "TEST( Feature, Covered ) {}\n",
        )
        head_sha = self.fixture.commit( "add unit test" )
        lcov = self.write_lcov( "" )
        report, markdown = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha ),
        )

        patch = report["patch"]
        self.assertEqual( patch["changed_files"], 0 )
        self.assertEqual( patch["changed_lines"], 0 )
        self.assertEqual( patch["excluded_changed_files"], 1 )
        self.assertEqual( patch["excluded_changed_lines"], 1 )
        self.assertIn( "test/example/support files", markdown )

    def test_markdown_escapes_paths_and_function_names( self ) -> None:
        path = "src/coreComponents/evil|<x>.cpp"
        base_sha, head_sha, source = self.make_added_source( path )
        lcov = self.write_lcov(
            f"SF:{source}\n"
            "FN:1,bad|<T>\nFNDA:0,bad|<T>\nDA:1,0\nend_of_record\n"
        )
        _, markdown = compare_pr_coverage.analyze(
            self.fixture.repository,
            base_sha,
            head_sha,
            lcov,
            self.candidate_summary( head_sha ),
        )

        self.assertIn( "evil&#124;&lt;x&gt;.cpp", markdown )
        self.assertIn( "bad&#124;&lt;T&gt;", markdown )
        self.assertNotIn( "evil|<x>.cpp", markdown )
        self.assertNotIn( "bad|<T>", markdown )

    def test_cli_writes_compact_outputs_and_rejects_short_sha( self ) -> None:
        base_sha, head_sha, source = self.make_added_source()
        lcov = self.write_lcov(
            f"SF:{source}\nDA:1,1\nend_of_record\n"
        )
        candidate = self.candidate_summary( head_sha )
        output_json = pathlib.Path( self.fixture.temporary.name ) / "diff.json"
        output_markdown = pathlib.Path( self.fixture.temporary.name ) / "diff.md"
        result = subprocess.run(
            [
                sys.executable,
                str( SCRIPT ),
                "--repository",
                str( self.fixture.repository ),
                "--base-sha",
                base_sha,
                "--head-sha",
                head_sha,
                "--native-info",
                str( lcov ),
                "--candidate-summary",
                str( candidate ),
                "--output-json",
                str( output_json ),
                "--output-markdown",
                str( output_markdown ),
            ],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual( result.returncode, 0, result.stderr )
        self.assertEqual( json.loads( output_json.read_text() )["schema_version"], 1 )
        self.assertIn( "Pull-request coverage", output_markdown.read_text() )
        self.assertLess( output_json.stat().st_size, 20_000 )
        self.assertEqual( stat.S_IMODE( output_json.stat().st_mode ), 0o644 )
        self.assertEqual( stat.S_IMODE( output_markdown.stat().st_mode ), 0o644 )

        invalid = subprocess.run(
            [
                sys.executable,
                str( SCRIPT ),
                "--repository",
                str( self.fixture.repository ),
                "--base-sha",
                base_sha[:12],
                "--head-sha",
                head_sha,
                "--native-info",
                str( lcov ),
                "--candidate-summary",
                str( candidate ),
            ],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual( invalid.returncode, 2 )
        self.assertIn( "full 40-character Git SHA", invalid.stderr )


if __name__ == "__main__":
    unittest.main()
