#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-2.1-only
#
# Merge Clang source-coverage profiles and report two production branch metrics:
#
#   canonical   LLVM's CoverageReport summary for the compiled configuration.
#   native      Every LCOV source-branch outcome emitted by llvm-cov. Template
#               and inline-function instantiations can appear more than once.
#
# All metrics include source files below src/coreComponents, excluding
# test/support directory names listed in production_exclude_regex below.

set -euo pipefail

# Coverage mapping never needs distribution debug packages. Disabling the
# host's optional debuginfod setting keeps reports local and reproducible.
export DEBUGINFOD_URLS=

usage()
{
  cat <<'EOF'
Usage: scripts/llvm_source_branch_coverage.sh BUILD_DIR PROFILE_DIR [OUTPUT_DIR]

BUILD_DIR    A CPU shared-library build configured with
             GEOS_ENABLE_LLVM_SOURCE_COVERAGE=ON.
PROFILE_DIR  Directory containing the *.profraw files from the test run.
OUTPUT_DIR   Report directory (default: BUILD_DIR/llvm-source-coverage).

Environment overrides:
  LLVM_COV       llvm-cov executable (default: compiler-matched llvm-cov)
  LLVM_PROFDATA  llvm-profdata executable (default: compiler-matched llvm-profdata)
  LLVM_READELF   llvm-readelf executable (default: compiler-matched llvm-readelf)
  GEOS_SOURCE_DIR
                 coreComponents source directory
                 (default: PROJECT_ROOT/src/coreComponents)

CI provenance (required):
  GEOS_CI_CONTAINER_IMAGE
  GEOS_CI_CONTAINER_IMAGE_ID
  GEOS_CI_CONTAINER_IMAGE_DIGESTS
  GEOS_CI_LLVM_PACKAGE_VERSIONS
                 Container label, immutable image ID, and JSON array of
                 immutable repository digests plus the exact runtime-installed
                 LLVM package revisions supplied by the CI launcher.

Typical profile collection:
  mkdir -p BUILD_DIR/profiles
  LLVM_PROFILE_FILE='BUILD_DIR/profiles/%p-%m.profraw' \
    ctest --test-dir BUILD_DIR --parallel 12 --output-on-failure

The script writes merged.profdata, llvm-summary.json, coverage-summary.json,
native.info, export logs, and branch-summary.txt. The JSON and .info files make
every numerator and denominator auditable with LLVM and LCOV tooling.
coverage-summary.json is the stable contract consumed by CI.
EOF
}

if (( $# < 2 || $# > 3 )); then
  usage >&2
  exit 2
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(dirname -- "${script_dir}")"
build_dir="$(cd -- "$1" && pwd -P)"
profile_dir="$(cd -- "$2" && pwd -P)"
output_arg="${3:-${build_dir}/llvm-source-coverage}"
mkdir -p -- "${output_arg}"
output_dir="$(cd -- "${output_arg}" && pwd -P)"
source_dir="${GEOS_SOURCE_DIR:-${project_root}/src/coreComponents}"
source_dir="$(cd -- "${source_dir}" && pwd -P)"
# A reused report directory must never mix fresh summaries with stale evidence
# from an earlier invocation that failed before exporting every format.
rm -f -- \
  "${output_dir}/branch-summary.txt" \
  "${output_dir}/branch-summary.txt.tmp" \
  "${output_dir}/coverage-objects.txt" \
  "${output_dir}/coverage-summary.json" \
  "${output_dir}/coverage-summary.json.tmp" \
  "${output_dir}/definition-export.log" \
  "${output_dir}/definition.info" \
  "${output_dir}/llvm-summary.json" \
  "${output_dir}/mapping-integrity.log" \
  "${output_dir}/merged.profdata" \
  "${output_dir}/native-export.log" \
  "${output_dir}/native.info" \
  "${output_dir}/profraw-inputs.txt" \
  "${output_dir}/summary-export.log"

cache_file="${build_dir}/CMakeCache.txt"
if [[ ! -f "${cache_file}" ]]; then
  echo "CMake cache not found: ${cache_file}" >&2
  exit 1
fi
if ! grep -q '^GEOS_ENABLE_LLVM_SOURCE_COVERAGE:BOOL=ON$' "${cache_file}"; then
  echo "Build does not have GEOS_ENABLE_LLVM_SOURCE_COVERAGE=ON" >&2
  exit 1
fi
if ! grep -q '^ENABLE_COVERAGE:BOOL=OFF$' "${cache_file}"; then
  echo "Build must have BLT ENABLE_COVERAGE=OFF" >&2
  exit 1
fi

c_compiler="$(sed -n 's/^CMAKE_C_COMPILER:[^=]*=//p' "${cache_file}" | head -n 1)"
cxx_compiler="$(sed -n 's/^CMAKE_CXX_COMPILER:[^=]*=//p' "${cache_file}" | head -n 1)"
for compiler in "${c_compiler}" "${cxx_compiler}"; do
  if [[ ! -x "${compiler}" ]]; then
    echo "Compiler from CMake cache is not executable: ${compiler}" >&2
    exit 1
  fi
done

c_compiler_major="$("${c_compiler}" --version | sed -nE '1s/.*version[[:space:]]+([0-9]+).*/\1/p')"
cxx_compiler_major="$("${cxx_compiler}" --version | sed -nE '1s/.*version[[:space:]]+([0-9]+).*/\1/p')"
if [[ -z "${c_compiler_major}" || -z "${cxx_compiler_major}" ||
      "${c_compiler_major}" != "${cxx_compiler_major}" ]]; then
  echo "Clang C/C++ compiler major versions are missing or inconsistent." >&2
  exit 1
fi
compiler_major="${cxx_compiler_major}"

resolve_llvm_tool()
{
  local environment_name="$1"
  local base_name="$2"
  local override="${!environment_name:-}"

  if [[ -n "${override}" ]]; then
    command -v -- "${override}"
  elif command -v -- "${base_name}-${compiler_major}" >/dev/null 2>&1; then
    command -v -- "${base_name}-${compiler_major}"
  else
    command -v -- "${base_name}"
  fi
}

llvm_cov="$(resolve_llvm_tool LLVM_COV llvm-cov || true)"
llvm_profdata="$(resolve_llvm_tool LLVM_PROFDATA llvm-profdata || true)"
llvm_readelf="$(resolve_llvm_tool LLVM_READELF llvm-readelf || true)"
python="${PYTHON:-python3}"

for tool in "${llvm_cov}" "${llvm_profdata}" "${llvm_readelf}" "${python}"; do
  if ! command -v -- "${tool}" >/dev/null 2>&1; then
    echo "Required tool not found: ${tool}" >&2
    exit 2
  fi
done

for tool in "${llvm_cov}" "${llvm_profdata}" "${llvm_readelf}"; do
  tool_major="$("${tool}" --version | sed -nE '1s/.*version[[:space:]]+([0-9]+).*/\1/p')"
  if [[ "${tool_major}" != "${compiler_major}" ]]; then
    echo "LLVM tool/compiler major-version mismatch: ${tool} is ${tool_major}, compiler is ${compiler_major}" >&2
    exit 2
  fi
done

profile_list="${output_dir}/profraw-inputs.txt"
find "${profile_dir}" -type f -name '*.profraw' -print | LC_ALL=C sort > "${profile_list}"
if [[ ! -s "${profile_list}" ]]; then
  echo "No .profraw files found below ${profile_dir}" >&2
  exit 1
fi

profdata="${output_dir}/merged.profdata"
"${llvm_profdata}" merge --input-files="${profile_list}" --output="${profdata}"

# The denominator is the compiled CPU product: geosx and its component shared
# libraries. Test executables are profile producers, not shipped artifacts;
# including their private header instantiations would make the production
# denominator depend on the shape of the test suite.
declare -a candidate_objects=( "${build_dir}/bin/geosx" )
if [[ ! -x "${candidate_objects[0]}" || ! -d "${build_dir}/lib" ]]; then
  echo "geosx executable or build library directory not found below ${build_dir}" >&2
  exit 1
fi
while IFS= read -r -d '' candidate; do
  case "$(basename -- "${candidate}")" in
    libgtest*|libtestingUtilities.*) continue ;;
  esac
  candidate_objects+=( "${candidate}" )
done < <(
  find "${build_dir}/lib" -maxdepth 1 -type f \
    \( -name '*.so' -o -name '*.so.*' -o -name '*.dylib' -o -name '*.dll' \) \
    -print0 | LC_ALL=C sort -z
)

object_list="${output_dir}/coverage-objects.txt"
: > "${object_list}"
for candidate in "${candidate_objects[@]}"; do
  if "${llvm_readelf}" --sections "${candidate}" 2>/dev/null | grep '__llvm_covmap' >/dev/null; then
    printf '%s\n' "${candidate}" >> "${object_list}"
  fi
done

if [[ ! -s "${object_list}" ]]; then
  echo "No objects with Clang coverage mappings found below ${build_dir}" >&2
  exit 1
fi
if ! grep -Eq '/geosx(\.exe)?$' "${object_list}"; then
  echo "Instrumented geosx executable not found below ${build_dir}" >&2
  exit 1
fi
if ! grep -Eq '/libmainInterface\.(so(\..*)?|dylib|dll)$' "${object_list}"; then
  echo "Instrumented mainInterface shared library not found below ${build_dir}" >&2
  exit 1
fi

mapfile -t coverage_objects < "${object_list}"
oldest_profile=""
while IFS= read -r profile; do
  if [[ -z "${oldest_profile}" || "${oldest_profile}" -nt "${profile}" ]]; then
    oldest_profile="${profile}"
  fi
done < "${profile_list}"
for object in "${coverage_objects[@]}"; do
  if [[ "${object}" -nt "${oldest_profile}" ]]; then
    echo "Coverage object is newer than at least one input profile: ${object}" >&2
    exit 1
  fi
done

primary_object="${coverage_objects[0]}"
declare -a object_args=()
for object in "${coverage_objects[@]:1}"; do
  object_args+=( -object "${object}" )
done

production_exclude_regex='/(unitTests|unitTestUtilities|integrationTests|tests|examples|benchmarks|testingUtilities)/'
production_exclude_regex+='|/codingUtilities/UnitTestUtilities[.]hpp$'

native_info="${output_dir}/native.info"
llvm_summary="${output_dir}/llvm-summary.json"
common_export_args=(
  -instr-profile "${profdata}"
  --sources "${source_dir}"
  -ignore-filename-regex "${production_exclude_regex}"
  -num-threads=1
)

integrity_log="${output_dir}/mapping-integrity.log"
if ! "${llvm_cov}" report "${primary_object}" "${object_args[@]}" \
    "${common_export_args[@]}" -dump > /dev/null 2> "${integrity_log}"; then
  cat "${integrity_log}" >&2
  exit 1
fi

mismatch_count="$(awk \
  '/warning: [0-9]+ functions? have mismatched data/ { total += $2 } END { print total + 0 }' \
  "${integrity_log}")"
mismatch_details="$(grep -c '^hash-mismatch:' "${integrity_log}" || true)"
nonzero_mismatches="$(grep '^hash-mismatch:' "${integrity_log}" | grep -vc 'hash = 0x0$' || true)"
if (( mismatch_count != mismatch_details || nonzero_mismatches != 0 )); then
  cat "${integrity_log}" >&2
  echo "Refusing to report coverage after a nonzero or unaudited function-hash mismatch." >&2
  exit 1
fi

integrity_warning_regex='hash mismatch|profile data.*out of date|object.*newer than profile'
integrity_warning_regex+='|no coverage data found|failed to load coverage'
if grep -Eiq "${integrity_warning_regex}" "${integrity_log}"; then
  cat "${integrity_log}" >&2
  echo "Refusing to report coverage after an LLVM profile-integrity warning." >&2
  exit 1
fi

run_export()
{
  local output_file="$1"
  local log_file="$2"
  shift 2

  if ! "${llvm_cov}" export "${primary_object}" "${object_args[@]}" \
      "${common_export_args[@]}" "$@" > "${output_file}" 2> "${log_file}"; then
    cat "${log_file}" >&2
    return 1
  fi

  local export_mismatch_count
  export_mismatch_count="$(awk \
    '/warning: [0-9]+ functions? have mismatched data/ { total += $2 } END { print total + 0 }' \
    "${log_file}")"
  if (( export_mismatch_count != mismatch_count )); then
    cat "${log_file}" >&2
    echo "LLVM export mismatch count differs from the audited mapping count." >&2
    return 1
  fi

  if grep -Eiq "${integrity_warning_regex}" "${log_file}"; then
    cat "${log_file}" >&2
    echo "Refusing to report coverage after an LLVM profile-integrity warning." >&2
    return 1
  fi
}

native_instantiation_args=()
llvm_cov_export_help="$("${llvm_cov}" export --help 2>&1)"
if grep -q -- '--unify-instantiations' <<< "${llvm_cov_export_help}"; then
  native_instantiation_args=( -unify-instantiations=false )
fi
run_export "${native_info}" "${output_dir}/native-export.log" \
  -format=lcov "${native_instantiation_args[@]}"
run_export "${llvm_summary}" "${output_dir}/summary-export.log" \
  -format=text -summary-only

summarize_branches()
{
  # Count the BRDA outcome records directly. This deliberately does not trust
  # BRF/BRH summaries because older llvm-cov versions can summarize template
  # instances differently from the records emitted in the same LCOV file.
  awk -F, '
    /^BRDA:/ {
      total += 1
      if ($4 != "-" && ($4 + 0) > 0) {
        covered += 1
      }
    }
    END {
      percent = total ? 100.0 * covered / total : 0.0
      printf "%d %d %.6f\n", covered, total, percent
    }
  ' "$1"
}

read -r native_covered native_total native_percent < <(summarize_branches "${native_info}")
read -r canonical_covered canonical_total canonical_percent < <(
  "${python}" - "${llvm_summary}" <<'PY'
import json
import sys

with open( sys.argv[1], encoding="utf-8" ) as stream:
  document = json.load( stream )

branches = document["data"][0]["totals"]["branches"]
covered = int( branches["covered"] )
total = int( branches["count"] )
percent = 100.0 * covered / total if total else 0.0
print( covered, total, f"{percent:.6f}" )
PY
)
profile_count="$(wc -l < "${profile_list}")"
object_count="${#coverage_objects[@]}"
zero_hash_mappings="${mismatch_details}"

if (( canonical_total == 0 || native_total == 0 )); then
  echo "LLVM reported a zero production branch denominator." >&2
  exit 1
fi

summary="${output_dir}/branch-summary.txt"
summary_tmp="${summary}.tmp"
{
  printf 'scope: %s\n' "${source_dir}"
  printf 'excluded directory regex: %s\n' "${production_exclude_regex}"
  printf 'profiles: %d\n' "${profile_count}"
  printf 'coverage objects: %d\n' "${object_count}"
  printf 'zero-hash mappings without profile counters: %d\n' "${zero_hash_mappings}"
  printf 'canonical LLVM branches: %d/%d (%.6f%%)\n' \
    "${canonical_covered}" "${canonical_total}" "${canonical_percent}"
  printf 'native LCOV branch outcomes: %d/%d (%.6f%%)\n' \
    "${native_covered}" "${native_total}" "${native_percent}"
} | tee "${summary_tmp}"
mv -- "${summary_tmp}" "${summary}"

coverage_summary="${output_dir}/coverage-summary.json"
coverage_summary_tmp="${coverage_summary}.tmp"
"${python}" - "${llvm_summary}" "${coverage_summary_tmp}" \
  "${project_root}" "${source_dir}" "${production_exclude_regex}" \
  "${compiler_major}" "${profile_count}" "${object_count}" \
  "${zero_hash_mappings}" \
  "${native_covered}" "${native_total}" \
  "${cache_file}" "${c_compiler}" "${cxx_compiler}" "${llvm_cov}" <<'PY'
import hashlib
import json
import os
import re
import subprocess
import sys

(
  llvm_summary_path,
  output_path,
  project_root,
  source_dir,
  excluded_regex,
  compiler_major,
  profile_count,
  object_count,
  zero_hash_mappings,
  native_covered,
  native_total,
  cache_file,
  c_compiler,
  cxx_compiler,
  llvm_cov,
) = sys.argv[1:]

COVERAGE_CONTRACT_ID = "geos-llvm-source-coverage-v1"
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


def command_output( description, command ):
  try:
    result = subprocess.run(
      command,
      check=True,
      stdout=subprocess.PIPE,
      stderr=subprocess.STDOUT,
      text=True,
    )
  except ( OSError, subprocess.CalledProcessError ) as error:
    output = getattr( error, "stdout", "" ) or ""
    raise ValueError(
      f"cannot determine {description}: {output.strip() or error}"
    ) from error
  output = result.stdout.rstrip( "\r\n" )
  if not output:
    raise ValueError( f"cannot determine {description}: empty output" )
  return output


def normalize_cmake_bool( value ):
  value = value.strip()
  upper = value.upper()
  if upper in ( "", "0", "FALSE", "IGNORE", "N", "NO", "NOTFOUND", "OFF" ):
    return False
  if upper.endswith( "-NOTFOUND" ):
    return False
  return True


def normalized_build_config( path ):
  entries = {}
  cache_entry = re.compile( r"^([^#/:=]+):([^=]+)=(.*)$" )
  with open( path, encoding="utf-8" ) as stream:
    for raw_line in stream:
      match = cache_entry.match( raw_line.rstrip( "\r\n" ) )
      if match is None:
        continue
      name, entry_type, value = match.groups()
      if entry_type == "INTERNAL":
        continue
      if name in entries:
        raise ValueError( f"duplicate CMake cache entry: {name}" )
      entries[name] = ( entry_type, value.strip() )

  if "CMAKE_BUILD_TYPE" not in entries:
    raise ValueError( "CMake cache does not define CMAKE_BUILD_TYPE" )
  build_type = entries["CMAKE_BUILD_TYPE"][1]
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
  fixed_names = {
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

  config = {}
  for name, ( entry_type, value ) in entries.items():
    relevant = (
      name in fixed_names
      or name in active_flag_names
      or name.startswith( "ENABLE_" )
      or name.startswith( "GEOS_BUILD_" )
      or name.startswith( "GEOS_ENABLE_" )
      or name.startswith( "LVARRAY_" )
      or name.startswith( "RAJA_ENABLE_" )
    )
    if not relevant:
      continue
    config[name] = normalize_cmake_bool( value ) if entry_type == "BOOL" else value

  required = {
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
  missing = sorted( required - config.keys() )
  if missing:
    raise ValueError(
      "CMake coverage configuration is incomplete: " + ", ".join( missing )
    )
  if config["ENABLE_COVERAGE"] is not False:
    raise ValueError( "ENABLE_COVERAGE must be OFF in normalized build config" )
  if config["GEOS_ENABLE_LLVM_SOURCE_COVERAGE"] is not True:
    raise ValueError(
      "GEOS_ENABLE_LLVM_SOURCE_COVERAGE must be ON in normalized build config"
    )
  if config["BUILD_SHARED_LIBS"] is not True or \
      config["GEOS_BUILD_SHARED_LIBS"] is not True:
    raise ValueError( "coverage object selection requires shared libraries" )
  return dict( sorted( config.items() ) )


def container_provenance():
  environment_names = (
    "GEOS_CI_CONTAINER_IMAGE",
    "GEOS_CI_CONTAINER_IMAGE_ID",
    "GEOS_CI_CONTAINER_IMAGE_DIGESTS",
  )
  missing = [ name for name in environment_names if not os.environ.get( name ) ]
  if missing:
    raise ValueError(
      "missing CI container provenance: " + ", ".join( missing )
    )

  image = os.environ["GEOS_CI_CONTAINER_IMAGE"]
  image_id = os.environ["GEOS_CI_CONTAINER_IMAGE_ID"]
  try:
    image_digests = json.loads( os.environ["GEOS_CI_CONTAINER_IMAGE_DIGESTS"] )
  except json.JSONDecodeError as error:
    raise ValueError( "GEOS_CI_CONTAINER_IMAGE_DIGESTS is not valid JSON" ) from error
  if (
    not isinstance( image_digests, list )
    or not image_digests
    or any( not isinstance( digest, str ) or not digest for digest in image_digests )
  ):
    raise ValueError(
      "GEOS_CI_CONTAINER_IMAGE_DIGESTS must be a nonempty JSON string array"
    )
  image_digests = sorted( set( image_digests ) )
  if re.fullmatch( r"sha256:[0-9a-f]{64}", image_id ) is None:
    raise ValueError( "GEOS_CI_CONTAINER_IMAGE_ID must be an immutable SHA256 ID" )
  for digest in image_digests:
    if re.fullmatch( r"[^\s@]+@sha256:[0-9a-f]{64}", digest ) is None:
      raise ValueError( "container repository digests must be immutable SHA256 digests" )
  if any( character in image for character in "\r\n\0" ):
    raise ValueError( "GEOS_CI_CONTAINER_IMAGE contains a control character" )
  return {
    "image": image,
    "image_id": image_id,
    "image_digests": image_digests,
  }


git_command = [
  "git",
  "-c",
  f"safe.directory={project_root}",
  "-C",
  project_root,
  "rev-parse",
  "--verify",
]
commit_sha = command_output( "Git HEAD commit", git_command + [ "HEAD^{commit}" ] )
tree_sha = command_output( "Git HEAD tree", git_command + [ "HEAD^{tree}" ] )
git_object_id = re.compile( r"(?:[0-9a-f]{40}|[0-9a-f]{64})" )
if git_object_id.fullmatch( commit_sha ) is None or \
    git_object_id.fullmatch( tree_sha ) is None:
  raise ValueError( "Git commit/tree provenance is not a canonical object ID" )

c_compiler_version = command_output( "C compiler version", [ c_compiler, "--version" ] )
cxx_compiler_version = command_output(
  "C++ compiler version", [ cxx_compiler, "--version" ]
)
llvm_cov_version = command_output( "llvm-cov version", [ llvm_cov, "--version" ] )
c_compiler_target = command_output( "C compiler target", [ c_compiler, "-dumpmachine" ] )
cxx_compiler_target = command_output(
  "C++ compiler target", [ cxx_compiler, "-dumpmachine" ]
)
if c_compiler_target != cxx_compiler_target:
  raise ValueError( "C and C++ compiler targets differ" )

container = container_provenance()
toolchain = {
  "c_compiler_version": c_compiler_version,
  "cxx_compiler_version": cxx_compiler_version,
  "llvm_cov_version": llvm_cov_version,
  "compiler_target": cxx_compiler_target,
}
package_versions = os.environ.get( "GEOS_CI_LLVM_PACKAGE_VERSIONS", "" ).splitlines()
if (
  not package_versions
  or package_versions != sorted( set( package_versions ) )
  or any( re.fullmatch( r"[A-Za-z0-9.+-]+=[^\s=]+", entry ) is None
          for entry in package_versions )
):
  raise ValueError(
    "GEOS_CI_LLVM_PACKAGE_VERSIONS must be a sorted package=version list"
  )
toolchain["llvm_package_versions"] = package_versions
build_config = normalized_build_config( cache_file )

with open( llvm_summary_path, encoding="utf-8" ) as stream:
  llvm_document = json.load( stream )

totals = llvm_document["data"][0]["totals"]

def metric( covered, total ):
  covered = int( covered )
  total = int( total )
  if covered < 0 or total < 0 or covered > total:
    raise ValueError( "LLVM emitted invalid coverage counts" )
  return {
    "covered": covered,
    "total": total,
    "not_covered": total - covered,
    "percent": round( 100.0 * covered / total if total else 0.0, 6 ),
  }

canonical = {}
for name in ( "regions", "functions", "lines", "branches" ):
  canonical[name] = metric( totals[name]["covered"], totals[name]["count"] )

source_root = os.path.realpath( source_dir )
branch_gaps = []
per_file_metrics = []
file_aggregates = {
  name: { "covered": 0, "total": 0 }
  for name in ( "regions", "functions", "lines", "branches" )
}
seen_paths = set()
for source in llvm_document["data"][0]["files"]:
  raw_filename = source["filename"]
  if not os.path.isabs( raw_filename ):
    raw_filename = os.path.join( project_root, raw_filename )
  filename = os.path.realpath( raw_filename )
  try:
    in_scope = os.path.commonpath( ( source_root, filename ) ) == source_root
  except ValueError:
    in_scope = False
  if not in_scope:
    continue
  relative_path = os.path.relpath( filename, project_root ).replace( os.sep, "/" )
  if relative_path in seen_paths:
    raise ValueError( f"duplicate LLVM per-file coverage path: {relative_path}" )
  seen_paths.add( relative_path )
  file_metrics = {}
  for name in ( "regions", "functions", "lines", "branches" ):
    llvm_metric = source["summary"][name]
    file_metrics[name] = metric( llvm_metric["covered"], llvm_metric["count"] )
    file_aggregates[name]["covered"] += file_metrics[name]["covered"]
    file_aggregates[name]["total"] += file_metrics[name]["total"]
  per_file_metrics.append( { "path": relative_path, "metrics": file_metrics } )

  branch_metric = file_metrics["branches"]
  if branch_metric["not_covered"] == 0:
    continue
  branch_gaps.append(
    {
      "path": relative_path,
      **branch_metric,
    }
  )
per_file_metrics.sort( key=lambda source: source["path"] )
branch_gaps.sort( key=lambda gap: ( -gap["not_covered"], gap["path"] ) )
for name, aggregate in file_aggregates.items():
  if ( aggregate["covered"], aggregate["total"] ) != (
    canonical[name]["covered"], canonical[name]["total"]
  ):
    raise ValueError( f"per-file {name} counts do not match canonical totals" )

scope = os.path.relpath( source_dir, project_root ).replace( os.sep, "/" )
contract_payload = {
  "summary_schema_version": 3,
  "contract_id": COVERAGE_CONTRACT_ID,
  "scope": scope,
  "excluded_regex": excluded_regex,
  "metric_semantics": METRIC_SEMANTICS,
  "container": {
    "image_id": container["image_id"],
    "image_digests": container["image_digests"],
  },
  "toolchain": toolchain,
  "build_config": build_config,
  "object_selection": OBJECT_SELECTION,
}
contract_bytes = json.dumps(
  contract_payload,
  allow_nan=False,
  ensure_ascii=True,
  separators=( ",", ":" ),
  sort_keys=True,
).encode( "utf-8" )
contract_fingerprint = hashlib.sha256( contract_bytes ).hexdigest()

document = {
  "schema_version": 3,
  "scope": scope,
  "excluded_regex": excluded_regex,
  "measurement": {
    "commit_sha": commit_sha,
    "tree_sha": tree_sha,
    "contract_id": COVERAGE_CONTRACT_ID,
    "contract_fingerprint": contract_fingerprint,
    "container": container,
    "toolchain": toolchain,
    "build_config": build_config,
    "metric_semantics": METRIC_SEMANTICS,
    "object_selection": OBJECT_SELECTION,
  },
  "tool": {
    "name": "llvm-cov",
    "major": int( compiler_major ),
  },
  "inputs": {
    "profiles": int( profile_count ),
    "coverage_objects": int( object_count ),
    "zero_hash_mappings": int( zero_hash_mappings ),
  },
  "metrics": canonical,
  "per_file_metrics": per_file_metrics,
  "top_branch_gaps": branch_gaps[:5],
  "supplemental": {
    "native_branch_outcomes": metric( native_covered, native_total ),
  },
}

with open( output_path, "w", encoding="utf-8" ) as stream:
  json.dump( document, stream, indent=2, sort_keys=True )
  stream.write( "\n" )
PY
mv -- "${coverage_summary_tmp}" "${coverage_summary}"
