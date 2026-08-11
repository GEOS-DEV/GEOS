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
  "${native_covered}" "${native_total}" <<'PY'
import json
import os
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
) = sys.argv[1:]

with open( llvm_summary_path, encoding="utf-8" ) as stream:
  llvm_document = json.load( stream )

totals = llvm_document["data"][0]["totals"]

def metric( covered, total ):
  covered = int( covered )
  total = int( total )
  return {
    "covered": covered,
    "total": total,
    "not_covered": total - covered,
    "percent": round( 100.0 * covered / total if total else 0.0, 6 ),
  }

canonical = {}
for name in ( "regions", "functions", "lines", "branches" ):
  canonical[name] = metric( totals[name]["covered"], totals[name]["count"] )

document = {
  "schema_version": 1,
  "scope": os.path.relpath( source_dir, project_root ),
  "excluded_regex": excluded_regex,
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
  "supplemental": {
    "native_branch_outcomes": metric( native_covered, native_total ),
  },
}

with open( output_path, "w", encoding="utf-8" ) as stream:
  json.dump( document, stream, indent=2, sort_keys=True )
  stream.write( "\n" )
PY
mv -- "${coverage_summary_tmp}" "${coverage_summary}"
