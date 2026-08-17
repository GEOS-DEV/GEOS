#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

"""Validate and enforce repository-owned LLVM coverage thresholds."""

import argparse
import hashlib
import html
import json
import math
import pathlib
import re
import sys


COVERAGE_SUMMARY_SCHEMA_VERSION = 3
COVERAGE_POLICY_SCHEMA_VERSION = 1
COVERAGE_CONTRACT_ID = "geos-llvm-source-coverage-v1"
CANONICAL_METRICS = ( "regions", "functions", "lines", "branches" )
DISPLAY_METRICS = ( "lines", "functions", "branches", "regions" )
SUPPLEMENTAL_METRICS = (
    ( "native_branch_outcomes", "Native branch outcomes" ),
)
if set( DISPLAY_METRICS ) != set( CANONICAL_METRICS ):
    raise RuntimeError( "display metrics must match canonical coverage metrics" )

METRIC_SEMANTICS = {
    "canonical_metrics": "llvm-cov-summary-instantiation-groups-v1",
    "native_branch_outcomes": "llvm-cov-lcov-emitted-brda-records-v1",
    "source_selection": "llvm-cov-sources-and-ignore-regex-v1",
}
OBJECT_SELECTION = {
    "contract_version": 1,
    "primary_object": "bin/geosx",
    "additional_object_directory": "lib",
    "additional_object_globs": [ "*.dll", "*.dylib", "*.so", "*.so.*" ],
    "excluded_library_basenames": [ "libgtest*", "libtestingUtilities.*" ],
    "required_object_basenames": [ "geosx", "libmainInterface" ],
    "coverage_mapping_section": "__llvm_covmap",
    "test_executables": "excluded",
}
MEASUREMENT_KEYS = {
    "commit_sha",
    "tree_sha",
    "contract_id",
    "contract_fingerprint",
    "container",
    "toolchain",
    "build_config",
    "metric_semantics",
    "object_selection",
}
CMAKE_FIXED_CONFIG_NAMES = {
    "BLT_CXX_STD",
    "BUILD_SHARED_LIBS",
    "CMAKE_BUILD_TYPE",
    "CMAKE_CUDA_ARCHITECTURES",
    "CMAKE_C_EXTENSIONS",
    "CMAKE_C_STANDARD",
    "CMAKE_C_STANDARD_REQUIRED",
    "CMAKE_CXX_EXTENSIONS",
    "CMAKE_CXX_STANDARD",
    "CMAKE_CXX_STANDARD_REQUIRED",
    "CMAKE_INTERPROCEDURAL_OPTIMIZATION",
    "CMAKE_POSITION_INDEPENDENT_CODE",
    "CMAKE_UNITY_BUILD",
    "GEOS_GLOBALINDEX_TYPE",
    "GEOS_GLOBALINDEX_TYPE_FLAG",
    "GEOS_LA_INTERFACE",
    "GEOS_LA_INTERFACE_HYPRE",
    "GEOS_LOCALINDEX_TYPE",
    "GEOS_LOCALINDEX_TYPE_FLAG",
}
CMAKE_REQUIRED_CONFIG_NAMES = {
    "BLT_CXX_STD",
    "BUILD_SHARED_LIBS",
    "CMAKE_BUILD_TYPE",
    "ENABLE_COVERAGE",
    "ENABLE_CUDA",
    "ENABLE_HIP",
    "ENABLE_HYPRE",
    "ENABLE_MPI",
    "ENABLE_OPENMP",
    "ENABLE_TRILINOS",
    "GEOS_BUILD_SHARED_LIBS",
    "GEOS_ENABLE_BOUNDS_CHECK",
    "GEOS_ENABLE_LLVM_SOURCE_COVERAGE",
    "GEOS_GLOBALINDEX_TYPE",
    "GEOS_LA_INTERFACE",
    "GEOS_LOCALINDEX_TYPE",
    "LVARRAY_BOUNDS_CHECK",
    "RAJA_ENABLE_CUDA",
    "RAJA_ENABLE_HIP",
    "RAJA_ENABLE_OPENMP",
}


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


def require_exact_keys( value: dict, expected: set[str], description: str ) -> None:
    actual = set( value )
    if actual != expected:
        missing = sorted( expected - actual )
        unexpected = sorted( actual - expected )
        details = []
        if missing:
            details.append( "missing " + ", ".join( missing ) )
        if unexpected:
            details.append( "unexpected " + ", ".join( unexpected ) )
        raise ValueError( f"{description} fields are invalid: {'; '.join( details )}" )


def require_string( value: object, description: str, allow_newlines=False ) -> str:
    if not isinstance( value, str ) or not value:
        raise ValueError( f"{description} must be a nonempty string" )
    forbidden = "\0" if allow_newlines else "\0\r\n"
    if any( character in value for character in forbidden ):
        raise ValueError( f"{description} contains a control character" )
    if value != value.rstrip():
        raise ValueError( f"{description} must not have trailing whitespace" )
    return value


def validate_repo_path( value: object, description: str ) -> str:
    path = require_string( value, description )
    normalized = pathlib.PurePosixPath( path )
    if (
        normalized.is_absolute()
        or normalized.as_posix() != path
        or ".." in normalized.parts
        or "\\" in path
    ):
        raise ValueError( f"{description} must be a normalized repository-relative path" )
    return path


def validate_scoped_path( value: object, scope: str, description: str ) -> str:
    path = validate_repo_path( value, description )
    scope_prefix = f"{scope}/"
    if not path.startswith( scope_prefix ) or len( path ) == len( scope_prefix ):
        raise ValueError( f"{description} is outside the scope" )
    return path


def validate_metric(
    metric: object,
    description: str,
    minimum_total: int = 1,
    allowed_extra_fields: set[str] | None = None,
) -> dict:
    if not isinstance( metric, dict ):
        raise ValueError( f"{description} must be an object" )
    expected_fields = { "covered", "total", "not_covered", "percent" }
    if allowed_extra_fields:
        expected_fields.update( allowed_extra_fields )
    require_exact_keys( metric, expected_fields, description )

    covered = require_integer( metric.get( "covered" ), f"{description}.covered" )
    total = require_integer(
        metric.get( "total" ), f"{description}.total", minimum_total
    )
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
    expected_percent = round( 100.0 * covered / total if total else 0.0, 6 )
    if abs( float( percent ) - expected_percent ) > 1.0e-6:
        raise ValueError( f"{description}.percent is inconsistent with its counts" )

    return {
        "covered": covered,
        "total": total,
        "not_covered": not_covered,
        "percent": expected_percent,
    }


def validate_container( raw_container: object ) -> dict:
    if not isinstance( raw_container, dict ):
        raise ValueError( "measurement.container must be an object" )
    require_exact_keys(
        raw_container,
        { "image", "image_id", "image_digests" },
        "measurement.container",
    )
    image = require_string( raw_container.get( "image" ), "measurement.container.image" )
    if any( character.isspace() for character in image ):
        raise ValueError( "measurement.container.image must be a container reference" )
    image_id = require_string(
        raw_container.get( "image_id" ), "measurement.container.image_id"
    )
    if re.fullmatch( r"sha256:[0-9a-f]{64}", image_id ) is None:
        raise ValueError( "measurement.container.image_id must be an immutable SHA256 ID" )

    raw_digests = raw_container.get( "image_digests" )
    if not isinstance( raw_digests, list ) or not raw_digests:
        raise ValueError(
            "measurement.container.image_digests must be a nonempty array"
        )
    digests = []
    for index, raw_digest in enumerate( raw_digests ):
        digest = require_string(
            raw_digest, f"measurement.container.image_digests[{index}]"
        )
        if re.fullmatch( r"[^\s@]+@sha256:[0-9a-f]{64}", digest ) is None:
            raise ValueError(
                "measurement.container.image_digests entries must be immutable SHA256 digests"
            )
        digests.append( digest )
    if digests != sorted( set( digests ) ):
        raise ValueError(
            "measurement.container.image_digests must be sorted and unique"
        )
    return { "image": image, "image_id": image_id, "image_digests": digests }


def tool_major_from_version( version: str, description: str ) -> int:
    match = re.search( r"\bversion\s+([0-9]+)(?:[.\s]|$)", version )
    if match is None:
        raise ValueError( f"{description} does not contain a version number" )
    return int( match.group( 1 ) )


def validate_toolchain( raw_toolchain: object, expected_major: int ) -> dict:
    if not isinstance( raw_toolchain, dict ):
        raise ValueError( "measurement.toolchain must be an object" )
    keys = {
        "c_compiler_version",
        "cxx_compiler_version",
        "llvm_cov_version",
        "llvm_package_versions",
        "compiler_target",
    }
    require_exact_keys( raw_toolchain, keys, "measurement.toolchain" )
    string_keys = keys - { "llvm_package_versions" }
    toolchain = {
        name: require_string(
            raw_toolchain.get( name ),
            f"measurement.toolchain.{name}",
            allow_newlines=name.endswith( "_version" ),
        )
        for name in string_keys
    }
    raw_packages = raw_toolchain.get( "llvm_package_versions" )
    if not isinstance( raw_packages, list ) or not raw_packages:
        raise ValueError(
            "measurement.toolchain.llvm_package_versions must be a nonempty array"
        )
    packages = []
    for index, raw_package in enumerate( raw_packages ):
        package = require_string(
            raw_package,
            f"measurement.toolchain.llvm_package_versions[{index}]",
        )
        if re.fullmatch( r"[A-Za-z0-9.+-]+=[^\s=]+", package ) is None:
            raise ValueError(
                "measurement.toolchain.llvm_package_versions entries must be package=version"
            )
        packages.append( package )
    if packages != sorted( set( packages ) ):
        raise ValueError(
            "measurement.toolchain.llvm_package_versions must be sorted and unique"
        )
    toolchain["llvm_package_versions"] = packages
    target = toolchain["compiler_target"]
    if any( character.isspace() for character in target ):
        raise ValueError( "measurement.toolchain.compiler_target is invalid" )
    for name in (
        "c_compiler_version",
        "cxx_compiler_version",
        "llvm_cov_version",
    ):
        if tool_major_from_version(
            toolchain[name], f"measurement.toolchain.{name}"
        ) != expected_major:
            raise ValueError(
                f"measurement.toolchain.{name} does not match tool.major"
            )
    return { name: toolchain[name] for name in sorted( toolchain ) }


def validate_build_config( raw_config: object ) -> dict:
    if not isinstance( raw_config, dict ) or not raw_config:
        raise ValueError( "measurement.build_config must be a nonempty object" )
    build_type = raw_config.get( "CMAKE_BUILD_TYPE" )
    if not isinstance( build_type, str ) or not build_type:
        raise ValueError(
            "measurement.build_config.CMAKE_BUILD_TYPE must be a nonempty string"
        )
    build_type_suffix = re.sub( r"[^A-Za-z0-9]", "_", build_type ).upper()
    active_flag_names = {
        "CMAKE_C_FLAGS",
        "CMAKE_CXX_FLAGS",
        "CMAKE_EXE_LINKER_FLAGS",
        "CMAKE_MODULE_LINKER_FLAGS",
        "CMAKE_SHARED_LINKER_FLAGS",
        f"CMAKE_C_FLAGS_{build_type_suffix}",
        f"CMAKE_CXX_FLAGS_{build_type_suffix}",
        f"CMAKE_EXE_LINKER_FLAGS_{build_type_suffix}",
        f"CMAKE_MODULE_LINKER_FLAGS_{build_type_suffix}",
        f"CMAKE_SHARED_LINKER_FLAGS_{build_type_suffix}",
    }

    config = {}
    for name, value in raw_config.items():
        if not isinstance( name, str ) or re.fullmatch( r"[A-Z][A-Z0-9_]*", name ) is None:
            raise ValueError( "measurement.build_config has an invalid CMake key" )
        relevant = (
            name in CMAKE_FIXED_CONFIG_NAMES
            or name in active_flag_names
            or name.startswith( "ENABLE_" )
            or name.startswith( "GEOS_BUILD_" )
            or name.startswith( "GEOS_ENABLE_" )
            or name.startswith( "LVARRAY_" )
            or name.startswith( "RAJA_ENABLE_" )
        )
        if not relevant:
            raise ValueError(
                f"measurement.build_config contains an irrelevant key: {name}"
            )
        if isinstance( value, bool ):
            config[name] = value
        elif isinstance( value, str ):
            if value != value.strip() or any(
                character in value for character in "\0\r\n"
            ):
                raise ValueError(
                    f"measurement.build_config.{name} is not normalized"
                )
            config[name] = value
        else:
            raise ValueError(
                f"measurement.build_config.{name} must be a string or boolean"
            )

    missing = sorted( CMAKE_REQUIRED_CONFIG_NAMES - config.keys() )
    if missing:
        raise ValueError(
            "measurement.build_config is missing required keys: "
            + ", ".join( missing )
        )
    if config["ENABLE_COVERAGE"] is not False:
        raise ValueError( "measurement.build_config.ENABLE_COVERAGE must be false" )
    if config["GEOS_ENABLE_LLVM_SOURCE_COVERAGE"] is not True:
        raise ValueError(
            "measurement.build_config.GEOS_ENABLE_LLVM_SOURCE_COVERAGE must be true"
        )
    if config["BUILD_SHARED_LIBS"] is not True or \
            config["GEOS_BUILD_SHARED_LIBS"] is not True:
        raise ValueError( "measurement.build_config must select shared libraries" )
    return dict( sorted( config.items() ) )


def contract_payload( scope: str, excluded_regex: str, measurement: dict ) -> dict:
    container = measurement["container"]
    return {
        "summary_schema_version": COVERAGE_SUMMARY_SCHEMA_VERSION,
        "contract_id": measurement["contract_id"],
        "scope": scope,
        "excluded_regex": excluded_regex,
        "metric_semantics": measurement["metric_semantics"],
        "container": {
            "image_id": container["image_id"],
            "image_digests": container["image_digests"],
        },
        "toolchain": measurement["toolchain"],
        "build_config": measurement["build_config"],
        "object_selection": measurement["object_selection"],
    }


def compute_contract_fingerprint(
    scope: str, excluded_regex: str, measurement: dict
) -> str:
    serialized = json.dumps(
        contract_payload( scope, excluded_regex, measurement ),
        allow_nan=False,
        ensure_ascii=True,
        separators=( ",", ":" ),
        sort_keys=True,
    ).encode( "utf-8" )
    return hashlib.sha256( serialized ).hexdigest()


def validate_measurement(
    raw_measurement: object, scope: str, excluded_regex: str, tool_major: int
) -> dict:
    if not isinstance( raw_measurement, dict ):
        raise ValueError( "coverage summary measurement must be an object" )
    require_exact_keys( raw_measurement, MEASUREMENT_KEYS, "measurement" )

    commit_sha = require_string(
        raw_measurement.get( "commit_sha" ), "measurement.commit_sha"
    )
    tree_sha = require_string(
        raw_measurement.get( "tree_sha" ), "measurement.tree_sha"
    )
    git_object_id = r"(?:[0-9a-f]{40}|[0-9a-f]{64})"
    if re.fullmatch( git_object_id, commit_sha ) is None:
        raise ValueError( "measurement.commit_sha must be a canonical Git object ID" )
    if re.fullmatch( git_object_id, tree_sha ) is None:
        raise ValueError( "measurement.tree_sha must be a canonical Git object ID" )

    contract_id = require_string(
        raw_measurement.get( "contract_id" ), "measurement.contract_id"
    )
    if contract_id != COVERAGE_CONTRACT_ID:
        raise ValueError(
            f"measurement.contract_id must be {COVERAGE_CONTRACT_ID!r}"
        )
    fingerprint = require_string(
        raw_measurement.get( "contract_fingerprint" ),
        "measurement.contract_fingerprint",
    )
    if re.fullmatch( r"[0-9a-f]{64}", fingerprint ) is None:
        raise ValueError(
            "measurement.contract_fingerprint must be a lowercase SHA256 digest"
        )

    container = validate_container( raw_measurement.get( "container" ) )
    toolchain = validate_toolchain(
        raw_measurement.get( "toolchain" ), tool_major
    )
    build_config = validate_build_config( raw_measurement.get( "build_config" ) )
    if raw_measurement.get( "metric_semantics" ) != METRIC_SEMANTICS:
        raise ValueError( "measurement.metric_semantics is not the v1 contract" )
    if raw_measurement.get( "object_selection" ) != OBJECT_SELECTION:
        raise ValueError( "measurement.object_selection is not the v1 contract" )

    measurement = {
        "commit_sha": commit_sha,
        "tree_sha": tree_sha,
        "contract_id": contract_id,
        "contract_fingerprint": fingerprint,
        "container": container,
        "toolchain": toolchain,
        "build_config": build_config,
        "metric_semantics": dict( METRIC_SEMANTICS ),
        "object_selection": dict( OBJECT_SELECTION ),
    }
    expected_fingerprint = compute_contract_fingerprint(
        scope, excluded_regex, measurement
    )
    if fingerprint != expected_fingerprint:
        raise ValueError(
            "measurement.contract_fingerprint does not match the coverage contract"
        )
    return measurement


def validate_per_file_metrics(
    raw_files: object, scope: str, canonical: dict
) -> list[dict]:
    if not isinstance( raw_files, list ) or not raw_files:
        raise ValueError( "per_file_metrics must be a nonempty array" )
    files = []
    seen_paths = set()
    aggregates = {
        name: { "covered": 0, "total": 0, "not_covered": 0 }
        for name in CANONICAL_METRICS
    }
    for index, raw_file in enumerate( raw_files ):
        if not isinstance( raw_file, dict ):
            raise ValueError( f"per_file_metrics[{index}] must be an object" )
        require_exact_keys(
            raw_file, { "path", "metrics" }, f"per_file_metrics[{index}]"
        )
        path = validate_scoped_path(
            raw_file.get( "path" ), scope, f"per_file_metrics[{index}].path"
        )
        if path in seen_paths:
            raise ValueError( "per_file_metrics paths must be unique" )
        seen_paths.add( path )
        raw_metrics = raw_file.get( "metrics" )
        if not isinstance( raw_metrics, dict ):
            raise ValueError( f"per_file_metrics[{index}].metrics must be an object" )
        require_exact_keys(
            raw_metrics,
            set( CANONICAL_METRICS ),
            f"per_file_metrics[{index}].metrics",
        )
        metrics = {
            name: validate_metric(
                raw_metrics[name],
                f"per_file_metrics[{index}].metrics.{name}",
                minimum_total=0,
            )
            for name in CANONICAL_METRICS
        }
        for name, metric in metrics.items():
            for count_name in ( "covered", "total", "not_covered" ):
                aggregates[name][count_name] += metric[count_name]
        files.append( { "path": path, "metrics": metrics } )

    if files != sorted( files, key=lambda source: source["path"] ):
        raise ValueError( "per_file_metrics must be sorted by path" )
    for name in CANONICAL_METRICS:
        for count_name in ( "covered", "total", "not_covered" ):
            if aggregates[name][count_name] != canonical[name][count_name]:
                raise ValueError(
                    f"per_file_metrics {name} {count_name} does not match canonical totals"
                )
    return files


def validate_branch_gaps( raw_gaps: object, scope: str ) -> list[dict]:
    if not isinstance( raw_gaps, list ) or len( raw_gaps ) > 5:
        raise ValueError( "top_branch_gaps must be an array with at most 5 entries" )

    gaps = []
    seen_paths = set()
    for index, raw_gap in enumerate( raw_gaps ):
        if not isinstance( raw_gap, dict ):
            raise ValueError( f"top_branch_gaps[{index}] must be an object" )
        path = validate_scoped_path(
            raw_gap.get( "path" ), scope, f"top_branch_gaps[{index}].path"
        )
        if path in seen_paths:
            raise ValueError( "top_branch_gaps paths must be normalized and unique" )
        seen_paths.add( path )

        gap = validate_metric(
            raw_gap,
            f"top_branch_gaps[{index}]",
            allowed_extra_fields={ "path" },
        )
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
    scope = validate_repo_path( document.get( "scope" ), "coverage summary scope" )

    require_exact_keys(
        document,
        {
            "schema_version",
            "scope",
            "excluded_regex",
            "measurement",
            "tool",
            "inputs",
            "metrics",
            "per_file_metrics",
            "top_branch_gaps",
            "supplemental",
        },
        "coverage summary",
    )
    excluded_regex = require_string(
        document.get( "excluded_regex" ), "coverage summary excluded_regex"
    )

    tool = document.get( "tool" )
    if not isinstance( tool, dict ) or tool.get( "name" ) != "llvm-cov":
        raise ValueError( "coverage summary tool.name must be llvm-cov" )
    require_exact_keys( tool, { "name", "major" }, "tool" )
    tool_major = require_integer( tool.get( "major" ), "tool.major", 1 )

    measurement = validate_measurement(
        document.get( "measurement" ), scope, excluded_regex, tool_major
    )

    raw_inputs = document.get( "inputs" )
    if not isinstance( raw_inputs, dict ):
        raise ValueError( "coverage summary inputs must be an object" )
    require_exact_keys(
        raw_inputs,
        { "profiles", "coverage_objects", "zero_hash_mappings" },
        "inputs",
    )
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
    require_exact_keys( raw_metrics, set( CANONICAL_METRICS ), "metrics" )
    metrics = {
        name: validate_metric( raw_metrics.get( name ), f"metrics.{name}" )
        for name in CANONICAL_METRICS
    }

    per_file_metrics = validate_per_file_metrics(
        document.get( "per_file_metrics" ), scope, metrics
    )

    raw_supplemental = document.get( "supplemental" )
    if not isinstance( raw_supplemental, dict ):
        raise ValueError( "coverage summary supplemental must be an object" )
    require_exact_keys(
        raw_supplemental,
        { name for name, _ in SUPPLEMENTAL_METRICS },
        "supplemental",
    )
    supplemental = {
        name: validate_metric( raw_supplemental.get( name ), f"supplemental.{name}" )
        for name, _ in SUPPLEMENTAL_METRICS
    }
    branch_gaps = validate_branch_gaps( document.get( "top_branch_gaps" ), scope )
    if sum( gap["total"] for gap in branch_gaps ) > metrics["branches"]["total"]:
        raise ValueError( "top_branch_gaps totals exceed canonical branch coverage" )
    expected_branch_gaps = []
    for source in per_file_metrics:
        branches = source["metrics"]["branches"]
        if branches["not_covered"]:
            expected_branch_gaps.append( { "path": source["path"], **branches } )
    expected_branch_gaps.sort(
        key=lambda gap: ( -gap["not_covered"], gap["path"] )
    )
    if branch_gaps != expected_branch_gaps[:5]:
        raise ValueError(
            "top_branch_gaps does not match canonical per-file branch metrics"
        )
    return {
        "scope": scope,
        "excluded_regex": excluded_regex,
        "measurement": measurement,
        "tool": { "name": "llvm-cov", "major": tool_major },
        "inputs": inputs,
        "metrics": metrics,
        "per_file_metrics": per_file_metrics,
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

    if args.markdown is not None:
        args.markdown.parent.mkdir( parents=True, exist_ok=True )
        args.markdown.write_text( markdown, encoding="utf-8" )
    else:
        print( markdown, end="" )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit( main() )
