#!/bin/bash
#
# Run the GEOS integrated tests locally.
#
# On a build configured with ENABLE_HYPREDRV=ON this additionally performs the
# hypredrive equivalence workflow: the suite is run once through hypredrive and once
# through the legacy hypre path (GEOS_HYPREDRV_FORCE_LEGACY=1) against the same
# baselines, and the per-solve linear-solver iteration counts of the two passes are
# compared so that solver-quality regressions cannot pass silently.
#
# Typical usage:
#
#   # Run everything (requires baselines, see --rebaseline below)
#   scripts/runIntegratedTests.sh -b <build-dir>
#
#   # Run only the hypre MGR / hypredrive equivalence tests
#   scripts/runIntegratedTests.sh -b <build-dir> -f _mgr
#
#   # First time on a machine without baselines: create them from this build.
#   # NOTE: baselines must be produced by a build WITHOUT hypredrive
#   # (ENABLE_HYPREDRV=OFF) or, equivalently, by a forced-legacy run.
#   scripts/runIntegratedTests.sh -b <build-dir> -f _mgr --rebaseline --force-legacy
#
#   # Then check the hypredrive path against those baselines
#   scripts/runIntegratedTests.sh -b <build-dir> -f _mgr --equivalence
#

set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GEOS_SRC_DIR="$( dirname "${SCRIPT_DIR}" )"

BUILD_DIR=""
FILTER=""
MACHINE="openmpi"
MAX_RANKS=16
MPIRUN=""
REBASELINE=false
FORCE_LEGACY=false
EQUIVALENCE=false
CLEAN=false
EXTRA_ARGS=()

usage()
{
    sed -n '2,/^$/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<EOF
Options:
  -b, --build-dir DIR   GEOS build directory (required)
  -f, --filter EXPR     only run tests whose name contains EXPR (e.g. _mgr)
  -n, --max-ranks N     max MPI ranks available to the suite (default 16); tests
                        needing more are skipped by ATS
  -m, --machine NAME    ATS machine (default openmpi; the built-in "generic"
                        machine is sequential-only and skips all parallel tests)
      --mpirun PATH     mpirun executable used by the openmpi machine
                        (default: the mpirun found on PATH)
      --rebaseline      rebaseline failed tests after the run
      --force-legacy    run with GEOS_HYPREDRV_FORCE_LEGACY=1 (legacy hypre path)
      --equivalence     run the suite twice (hypredrive, then forced legacy) and
                        compare linear-solver iteration counts between the passes
      --clean           remove previous test results before running
  -h, --help            show this help
Any remaining arguments are forwarded to geos_ats.sh.
EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--build-dir)  BUILD_DIR=$2; shift 2;;
        -f|--filter)     FILTER=$2; shift 2;;
        -n|--max-ranks)  MAX_RANKS=$2; shift 2;;
        -m|--machine)    MACHINE=$2; shift 2;;
        --mpirun)        MPIRUN=$2; shift 2;;
        --rebaseline)    REBASELINE=true; shift;;
        --force-legacy)  FORCE_LEGACY=true; shift;;
        --equivalence)   EQUIVALENCE=true; shift;;
        --clean)         CLEAN=true; shift;;
        -h|--help)       usage; exit 0;;
        *)               EXTRA_ARGS+=("$1"); shift;;
    esac
done

if [[ -z ${BUILD_DIR} ]]; then
    echo "ERROR: a build directory is required (-b/--build-dir)."
    usage
    exit 1
fi
BUILD_DIR=$( cd "${BUILD_DIR}" && pwd ) || exit 1

ATS_SCRIPT=${BUILD_DIR}/integratedTests/geos_ats.sh
if [[ ! -x ${ATS_SCRIPT} ]]; then
    echo "ERROR: ${ATS_SCRIPT} not found."
    echo "Run 'make ats_environment' in the build directory first."
    echo "If that fails with 'externally-managed-environment', reconfigure the build"
    echo "with -DPython3_EXECUTABLE=<path to a virtualenv python>."
    exit 1
fi

# The generic ATS machine is sequential-only and skips every test needing more than
# one rank, so select a machine that schedules MPI ranks and cap it at MAX_RANKS.
# Machine options must be passed through --ats (see ci_build_and_test_in_container.sh);
# passing them as plain arguments leaves the machine unconfigured and the run hangs.
# openmpi_mpirun is given as an absolute path so it does not depend on an install prefix.
if [[ -z ${MPIRUN} ]]; then
    MPIRUN=$( command -v mpirun )
    if [[ -z ${MPIRUN} ]]; then
        echo "ERROR: no mpirun found on PATH; pass --mpirun <path>."
        exit 1
    fi
fi
ATS_ARGS=(--machine "${MACHINE}")
if [[ ${MACHINE} == openmpi ]]; then
    ATS_ARGS+=(--ats "openmpi_mpirun=${MPIRUN}"
               --ats "openmpi_maxprocs=${MAX_RANKS}"
               --ats "openmpi_procspernode=${MAX_RANKS}")
fi
# ATS filters are evaluated against each test's options; geos_ats labels its steps
# "<testCaseName>_<step>_<stepType>", so filtering on the label selects whole tests.
[[ -n ${FILTER} ]] && ATS_ARGS+=(--ats "filter='${FILTER}' in label")
ATS_ARGS+=("${EXTRA_ARGS[@]}")

runSuite()  # $1 = label, $2 = "legacy" to force the legacy hypre path
{
    local label=$1 mode=$2
    echo
    echo "=================================================================="
    echo "Integrated tests: ${label}"
    echo "=================================================================="
    if [[ ${mode} == legacy ]]; then
        GEOS_HYPREDRV_FORCE_LEGACY=1 "${ATS_SCRIPT}" "${ATS_ARGS[@]}"
    else
        "${ATS_SCRIPT}" "${ATS_ARGS[@]}"
    fi
    return $?
}

harvestIterations()  # $1 = output json
{
    python3 "${SCRIPT_DIR}/compareLinearSolverIterations.py" harvest \
            "${BUILD_DIR}/integratedTests/TestResults/test_data" -o "$1"
}

if [[ ${CLEAN} == true ]]; then
    echo "Cleaning previous test results..."
    "${ATS_SCRIPT}" -a veryclean > /dev/null 2>&1
fi

STATUS=0

if [[ ${EQUIVALENCE} == true ]]; then
    ITER_DIR=${BUILD_DIR}/integratedTests
    runSuite "hypredrive path" || STATUS=$?
    harvestIterations "${ITER_DIR}/iterations_hypredrive.json"

    # Confirm the hypredrive path was actually taken (guards against a silently
    # disabled ENABLE_HYPREDRV or a silent fallback to the legacy solver).
    if grep -rql "hypredrive input" "${BUILD_DIR}/integratedTests/TestResults" 2>/dev/null; then
        echo "Confirmed: hypredrive input banner found in the test logs."
    else
        echo "ERROR: no 'hypredrive input' banner in any test log."
        echo "       The hypredrive path was not exercised (is ENABLE_HYPREDRV ON?)."
        STATUS=1
    fi
    if grep -rql "falling back to the legacy hypre solver" "${BUILD_DIR}/integratedTests/TestResults" 2>/dev/null; then
        echo "WARNING: at least one solver fell back from hypredrive to legacy hypre:"
        grep -rl "falling back to the legacy hypre solver" \
             "${BUILD_DIR}/integratedTests/TestResults" 2>/dev/null | head -5
    fi

    "${ATS_SCRIPT}" -a veryclean > /dev/null 2>&1
    runSuite "legacy hypre path" legacy || STATUS=$?
    harvestIterations "${ITER_DIR}/iterations_legacy.json"

    echo
    echo "=================================================================="
    echo "Linear-solver iteration parity"
    echo "=================================================================="
    python3 "${SCRIPT_DIR}/compareLinearSolverIterations.py" compare \
            "${ITER_DIR}/iterations_hypredrive.json" \
            "${ITER_DIR}/iterations_legacy.json" || STATUS=1
else
    if [[ ${FORCE_LEGACY} == true ]]; then
        runSuite "legacy hypre path" legacy || STATUS=$?
    else
        runSuite "default solver path" || STATUS=$?
    fi
fi

if [[ ${REBASELINE} == true ]]; then
    echo
    echo "Rebaselining failed tests..."
    # geos_ats prompts per test and the flag is not settable from the command line,
    # so answer the prompts on stdin.
    yes y | "${ATS_SCRIPT}" "${ATS_ARGS[@]}" -a rebaselinefailed
fi

echo
echo "Test results: ${BUILD_DIR}/integratedTests/TestResults"
exit ${STATUS}
