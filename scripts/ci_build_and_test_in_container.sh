#!/bin/bash
set -o pipefail

export PYTHONDONTWRITEBYTECODE=1

SCRIPT_NAME=$0
echo "Running CLI ${SCRIPT_NAME} $@"


# docs.docker.com/config/containers/resource_constraints
# Inside the container, tools like free report the host's available swap, not what's available inside the container.
# Don't rely on the output of free or similar tools to determine whether swap is present.
echo "running free -g"
free -g

# The or_die function run the passed command line and
# exits the program in case of non zero error code
function or_die () {
    "$@"
    local status=$?

    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

tempdir=""

function openssl_fips_provider_available () {
    local fips_module_path

    for fips_module_path in \
        /usr/lib/x86_64-linux-gnu/ossl-modules/fips.so \
        /usr/lib/aarch64-linux-gnu/ossl-modules/fips.so \
        /usr/lib64/ossl-modules/fips.so \
        /usr/lib/ossl-modules/fips.so \
        /usr/local/lib64/ossl-modules/fips.so \
        /usr/local/lib/ossl-modules/fips.so; do
        if [[ -f "${fips_module_path}" ]]; then
            return 0
        fi
    done

    return 1
}

function configure_openssl_for_non_fips_ubuntu_container () {
    local fips_enabled=""
    local openssl_conf_path=/tmp/geos-openssl-non-fips.cnf

    if [[ -r /proc/sys/crypto/fips_enabled ]]; then
        fips_enabled="$(tr -d '[:space:]' < /proc/sys/crypto/fips_enabled)"
    fi

    if [[ "${fips_enabled}" != "1" ]]; then
        return 0
    fi

    if [[ ! -r /etc/os-release ]] || ! grep -Eq '^ID="?ubuntu"?$' /etc/os-release; then
        return 0
    fi

    if [[ -e /etc/system-fips ]]; then
        return 0
    fi

    if openssl_fips_provider_available; then
        return 0
    fi

    export OPENSSL_FORCE_FIPS_MODE=0

    if [[ -n "${OPENSSL_CONF:-}" ]]; then
        echo "Host FIPS mode is visible in this Ubuntu container, but OPENSSL_CONF is already set to ${OPENSSL_CONF}; leaving it unchanged."
        echo "Using OPENSSL_FORCE_FIPS_MODE=0 because no OpenSSL FIPS provider module was found."
        return 0
    fi

    if ! cat > "${openssl_conf_path}" <<'EOF'
openssl_conf = openssl_init

[openssl_init]
providers = provider_sect
alg_section = algorithm_sect

[provider_sect]
default = default_sect

[default_sect]
activate = 1

[algorithm_sect]
default_properties = fips=no
EOF
    then
        echo "WARNING: unable to write ${openssl_conf_path}; leaving OpenSSL configuration unchanged."
        return 0
    fi

    export OPENSSL_CONF="${openssl_conf_path}"

    echo "Host FIPS mode is visible in this non-FIPS Ubuntu container, and no OpenSSL FIPS provider module was found."
    echo "Using OPENSSL_CONF=${OPENSSL_CONF} so OpenSSL initializes the default provider for CI tooling."
}

function print_crypto_diagnostics () {
    local fips_enabled=""

    if [[ -r /proc/sys/crypto/fips_enabled ]]; then
        fips_enabled="$(tr -d '[:space:]' < /proc/sys/crypto/fips_enabled)"
    fi

    if [[ "${fips_enabled}" != "1" && -z "${OPENSSL_CONF:-}" && -z "${OPENSSL_FORCE_FIPS_MODE:-}" ]]; then
        return 0
    fi

    echo "Crypto/FIPS summary:"
    echo "  host fips_enabled: ${fips_enabled:-unavailable}"
    echo "  container /etc/system-fips: $(if [[ -e /etc/system-fips ]]; then echo present; else echo absent; fi)"
    echo "  openssl fips provider module: $(if openssl_fips_provider_available; then echo present; else echo absent; fi)"
    echo "  OPENSSL_CONF: ${OPENSSL_CONF:-<unset>}"
    echo "  OPENSSL_FORCE_FIPS_MODE: ${OPENSSL_FORCE_FIPS_MODE:-<unset>}"

    if command -v openssl >/dev/null 2>&1; then
        openssl version 2>&1 | sed 's/^/  openssl version: /' || true
    fi
}

function bootstrap_pip () {
    local python_executable="$1"

    if [[ "${GEOS_SKIP_PIP_BOOTSTRAP:-false}" == "true" ]]; then
        echo "Skipping pip bootstrap because GEOS_SKIP_PIP_BOOTSTRAP=true."
        return 0
    fi

    echo "Updating pip"
    if ! "${python_executable}" -m pip install --upgrade pip setuptools wheel; then
        echo "::warning::pip bootstrap failed for ${python_executable}; continuing with the existing pip so the package install step can report the actionable failure."
        "${python_executable}" -m pip --version || true
        return 0
    fi

    "${python_executable}" -m pip cache purge || true
}

function usage () {
>&2 cat << EOF
Usage: $0
  --build-exe-only
      Request for the build of geos only.
  --cmake-build-type ...
      One of Debug, Release, RelWithDebInfo and MinSizeRel. Forwarded to CMAKE_BUILD_TYPE.
  --cmake-cuda-architectures ...
      Optional override for CMAKE_CUDA_ARCHITECTURES.
  --coverage-output-dir /path/to/output
      Directory where LLVM source-coverage inputs and the HTML report are written.
  --coverage-base-sha SHA
      Exact pull-request base commit used for advisory patch coverage.
  --coverage-head-sha SHA
      Exact commit being built and measured by the coverage job.
  --coverage-baseline-summary /path/to/coverage-summary.json
      Optional trusted summary for the exact base commit.
  --ctest-parallel-level N
      Number of tests ctest may run in parallel.
  --data-basename output.tar.gz
      If some data needs to be extracted from the build, the argument will define the tarball. Has to be a `tar.gz`.
  --geos-enable-bounds-check
      Either ON or OFF (default is ON). Build geos with bounds check.
  --enable-hypre
      One of ON or OFF (default is ON). Build geos with hypre.
  --enable-hypredrv
      One of ON or OFF (default is ON). Build geos with hypredrive.
      This flag overrides the TPL image host-config.
  --enable-hypre-device
      One of CPU, CUDA, or HIP (default is CPU). Build geos with hypre GPU support.
  --enable-trilinos
      One of ON or OFF (default is OFF). Build geos with trilinos.
  --exchange-dir /path/to/exchange
      Folder to share data with outside of the container.
  --host-config host-config/my_config.cmake
      The host-config. Path is relative to the root of the repository.
  --install-dir-basename GEOS-e42ffc1
      GEOS installation basename.
  --llvm-source-coverage
      Build with Clang source coverage, run the coverage smoke suite, and
      enforce the repository coverage policy.
  --makefile
      Use "Unix Makefiles" as build system generator.
  --ninja
      Use "Ninja" as build system generator.
  --no-install-schema
      Do not install the xsd schema.
  --no-run-unit-tests
      Do not run the unit tests (but they will be built).
  --nproc N
      Number of cores to use for the build.
  --use-native-architecture
      Build with compiler flags targeting the native runner CPU.
  --repository /path/to/repository
      Internal mountpoint where the geos repository will be available.
  --run-integrated-tests
      Run the integrated tests. Then bundle and send the results to the cloud.
  --use-sccache
      Enable sccache as compiler launcher.
  --sccache-credentials credentials.json
      Basename of the json credentials file to connect to the sccache cloud cache.
  --test-code-style
  --test-documentation
  -h | --help
EOF
exit 1
}

# First working in the root of the cloned repository.
# Then we'll move to the build dir.
or_die cd $(dirname $0)/..

args=$(or_die getopt -a -o h --long build-exe-only,cmake-build-type:,cmake-cuda-architectures:,coverage-base-sha:,coverage-baseline-summary:,coverage-head-sha:,coverage-output-dir:,ctest-parallel-level:,data-basename:,geos-enable-bounds-check:,enable-hypre:,enable-hypredrv:,enable-hypre-device:,enable-trilinos:,exchange-dir:,host-config:,install-dir-basename:,llvm-source-coverage,makefile,ninja,no-install-schema,no-run-unit-tests,nproc:,repository:,run-integrated-tests,sccache-credentials:,test-code-style,test-documentation,use-native-architecture,use-sccache,help -- "$@")

# Variables with default values
BUILD_EXE_ONLY=false
BUILD_GENERATOR=""
GEOS_INSTALL_SCHEMA=true
HOST_CONFIG="host-configs/environment.cmake"
ENABLE_HYPRE=ON
ENABLE_HYPREDRV=ON
ENABLE_HYPRE_DEVICE=CPU
GEOS_LA_INTERFACE=Hypre
RUN_UNIT_TESTS=true
RUN_INTEGRATED_TESTS=false
UPLOAD_TEST_BASELINES=false
TEST_CODE_STYLE=false
TEST_DOCUMENTATION=false
ENABLE_TRILINOS=OFF
LLVM_SOURCE_COVERAGE=false
COVERAGE_OUTPUT_DIR=""
COVERAGE_BASE_SHA=""
COVERAGE_HEAD_SHA=""
COVERAGE_BASELINE_SUMMARY=""
CTEST_PARALLEL_LEVEL_ARG=""
NPROC="$(nproc)"
GEOS_ENABLE_BOUNDS_CHECK=ON
SCCACHE_BIN=""
SCCACHE_CREDS=""
USE_SCCACHE=false
CMAKE_CUDA_ARCHITECTURES_ARGS=()
CMAKE_NATIVE_ARCHITECTURE_ARGS=()
ATS_CMAKE_ARGS=()

eval set -- ${args}
while :
do
  case $1 in
    --build-exe-only)
      BUILD_EXE_ONLY=true
      RUN_UNIT_TESTS=false
      shift;;
    --cmake-build-type)      CMAKE_BUILD_TYPE=$2;        shift 2;;
    --cmake-cuda-architectures)
      CMAKE_CUDA_ARCHITECTURES_ARGS+=("-DCMAKE_CUDA_ARCHITECTURES=$2")
      shift 2;;
    --coverage-base-sha)      COVERAGE_BASE_SHA=$2;       shift 2;;
    --coverage-baseline-summary) COVERAGE_BASELINE_SUMMARY=$2; shift 2;;
    --coverage-head-sha)      COVERAGE_HEAD_SHA=$2;       shift 2;;
    --coverage-output-dir)    COVERAGE_OUTPUT_DIR=$2;     shift 2;;
    --ninja)
        BUILD_GENERATOR=$1;
        shift;;
    --data-basename)
      DATA_BASENAME=$2
      DATA_BASENAME_WE=${DATA_BASENAME%%.*}
      DATA_BASENAME_EXT=${DATA_BASENAME#*.}
      if [[ ${DATA_BASENAME_EXT} != "tar.gz" ]] ; then
          echo "The script ${SCRIPT_NAME} can only pack data into a '.tar.gz' file."
          exit 1
      fi
      unset DATA_BASENAME DATA_BASENAME_EXT
      shift 2;;
    --geos-enable-bounds-check) GEOS_ENABLE_BOUNDS_CHECK=$2; shift 2;;
    --enable-hypre)          ENABLE_HYPRE=$2;            shift 2;;
    --enable-hypredrv)       ENABLE_HYPREDRV=$2;         shift 2;;
    --enable-hypre-device)   ENABLE_HYPRE_DEVICE=$2;     shift 2;;
    --enable-trilinos)       ENABLE_TRILINOS=$2;         shift 2;;
    --exchange-dir)          DATA_EXCHANGE_DIR=$2;       shift 2;;
    --host-config)           HOST_CONFIG=$2;             shift 2;;
    --install-dir-basename)  GEOS_DIR=${GEOSX_TPL_DIR}/../$2; shift 2;;
    --llvm-source-coverage)  LLVM_SOURCE_COVERAGE=true;  shift;;
    --makefile)              BUILD_GENERATOR="";         shift;;
    --no-install-schema)     GEOS_INSTALL_SCHEMA=false; shift;;
    --no-run-unit-tests)     RUN_UNIT_TESTS=false;       shift;;
    --nproc)                 NPROC=$2;                   shift 2;;
    --use-native-architecture)
      CMAKE_NATIVE_ARCHITECTURE_ARGS+=('-DCMAKE_C_FLAGS:STRING="-march=native -mtune=native"')
      CMAKE_NATIVE_ARCHITECTURE_ARGS+=('-DCMAKE_CXX_FLAGS:STRING="-march=native -mtune=native"')
      shift;;
    --repository)            GEOS_SRC_DIR=$2;            shift 2;;
    --run-integrated-tests)  RUN_INTEGRATED_TESTS=true;  shift;;
    --upload-test-baselines) UPLOAD_TEST_BASELINES=true; shift;;
    --ctest-parallel-level)
      CTEST_PARALLEL_LEVEL_ARG=$2
      shift 2;;
    --sccache-credentials)   SCCACHE_CREDS=$2; USE_SCCACHE=true; shift 2;;
    --use-sccache)           USE_SCCACHE=true;           shift;;
    --test-code-style)       TEST_CODE_STYLE=true;       shift;;
    --test-documentation)    TEST_DOCUMENTATION=true;    shift;;
    -h | --help)             usage;                      shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

if [[ -z "${GEOS_SRC_DIR}" ]]; then
  echo "Variable GEOS_SRC_DIR is either empty or not defined. Please define it using '--repository'."
  exit 1
fi


cleanup() {
  echo "Container cleanup..."
  if [[ -n "${tempdir:-}" ]]; then
    rm -rf "${tempdir}" || true
  fi
  rm -rf "${GEOS_SRC_DIR}/src/docs/sphinx/datastructure" || true

  if [[ -n "${HOST_UID:-}" && -n "${HOST_GID:-}" ]]; then
    chown -R "${HOST_UID}:${HOST_GID}" "${GEOS_SRC_DIR}" || true
  fi
}
trap cleanup EXIT



if [[ -z "${GEOS_DIR}" ]]; then
  echo "Installation folder undefined. Set to default value '/dev/null'. You can define it using '--install-dir-basename'."
  GEOS_DIR=/dev/null
fi

if [[ "${ENABLE_HYPRE}" = ON ]]; then
  GEOS_LA_INTERFACE=Hypre
else
  GEOS_LA_INTERFACE=Trilinos
fi

configure_openssl_for_non_fips_ubuntu_container
print_crypto_diagnostics

if [[ "${USE_SCCACHE}" == true ]]; then
  SCCACHE_BIN=${SCCACHE:-$(command -v sccache || true)}

  if [[ -z "${SCCACHE_BIN}" ]]; then
    echo "sccache was requested, but no sccache binary is available in the container."
    exit 1
  fi

  if [[ -n "${SCCACHE_CREDS:-}" ]]; then
    # The credential json file is available at the root of the geos repository.
    # We hereafter create the config file that points to it.
    # We use this file since it's managed by the 'google-github-actions/auth' actions.
    or_die mkdir -p ${HOME}/.config/sccache
    or_die cat <<EOT >> ${HOME}/.config/sccache/config
[cache.gcs]
rw_mode = "READ_WRITE"
cred_path = "${GEOS_SRC_DIR}/${SCCACHE_CREDS}"
bucket = "geos-dev"
key_prefix = "sccache"
EOT
  fi

  # Backend-specific credentials and endpoints are injected through the environment or generated config.
  SCCACHE_CMAKE_ARGS="-DCMAKE_C_COMPILER_LAUNCHER=${SCCACHE_BIN} -DCMAKE_CXX_COMPILER_LAUNCHER=${SCCACHE_BIN} -DCMAKE_CUDA_COMPILER_LAUNCHER=${SCCACHE_BIN}"

  if [[ -f /certs/ca-bundle.crt ]]; then
    export SSL_CERT_FILE=/certs/ca-bundle.crt
    export CURL_CA_BUNDLE=/certs/ca-bundle.crt
    export REQUESTS_CA_BUNDLE=/certs/ca-bundle.crt
  fi

  echo "sccache enabled: ${SCCACHE_BIN}, cache directory: ${SCCACHE_DIR:-default}"
fi

if [ -z "${NPROC}" ]; then
  NPROC=$(nproc)
  echo "NPROC unset, setting to ${NPROC}..."
fi
if ! [[ "${NPROC}" =~ ^[1-9][0-9]*$ ]]; then
  echo "NPROC must be a positive integer, got '${NPROC}'."
  exit 1
fi
echo "Using ${NPROC} cores."

if [[ -n "${CTEST_PARALLEL_LEVEL_ARG}" ]]; then
  export CTEST_PARALLEL_LEVEL="${CTEST_PARALLEL_LEVEL_ARG}"
  echo "Running ctest with CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL}."
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  echo "Running the integrated tests has been requested."

  # We install the python environment required by ATS to run the integrated tests.
  ATS_PYTHON_HOME=/tmp/run_integrated_tests_virtualenv
  if ! python3 -m venv --system-site-packages "${ATS_PYTHON_HOME}" ||
     ! "${ATS_PYTHON_HOME}/bin/python3" -c 'import setuptools, wheel' >/dev/null 2>&1; then
    echo "ATS virtualenv lacks required Python packaging modules; installing distro packages and recreating it."
    rm -rf "${ATS_PYTHON_HOME}"
    or_die apt-get update
    or_die apt-get install -y python3-dev python3-venv python3-setuptools python3-wheel
    or_die python3 -m venv --system-site-packages "${ATS_PYTHON_HOME}"
    or_die "${ATS_PYTHON_HOME}/bin/python3" -c 'import setuptools, wheel'
  fi

  bootstrap_pip "${ATS_PYTHON_HOME}/bin/python3"

  # Setup a temporary directory to hold tests
  tempdir=$(mktemp -d)
  echo "Setting up a temporary directory to hold tests and baselines: $tempdir"
  ATS_BASELINE_DIR=$tempdir/GEOS_integratedTests_baselines
  ATS_WORKING_DIR=$tempdir/GEOS_integratedTests_working

  export ATS_FILTER="np<=32"
  # Open MPI refuses root launches in CI containers unless both guard variables are set.
  # Keep this separate from openmpi_args because repeated geos-ats overrides replace values.
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
  ATS_ARGUMENTS="--machine openmpi --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_procspernode=${NPROC} --ats openmpi_maxprocs=${NPROC} --ats cutoff=45m"
  echo "Integrated test ATS arguments: ${ATS_ARGUMENTS}"
  ATS_CMAKE_ARGS=("-DATS_ARGUMENTS:STRING=\"${ATS_ARGUMENTS}\""
                  "-DPython3_ROOT_DIR=${ATS_PYTHON_HOME}"
                  "-DPython3_EXECUTABLE=${ATS_PYTHON_HOME}/bin/python3"
                  "-DATS_BASELINE_DIR=${ATS_BASELINE_DIR}"
                  "-DATS_WORKING_DIR=${ATS_WORKING_DIR}")
fi


if [[ "${LLVM_SOURCE_COVERAGE}" = true ]]; then
  if [[ -z "${COVERAGE_OUTPUT_DIR}" ]]; then
    echo "--llvm-source-coverage requires --coverage-output-dir." >&2
    exit 1
  fi
  if [[ ! "${COVERAGE_HEAD_SHA}" =~ ^[0-9a-f]{40}$ ||
        "$(git -c "safe.directory=${GEOS_SRC_DIR}" -C "${GEOS_SRC_DIR}" rev-parse HEAD)" != "${COVERAGE_HEAD_SHA}" ]]; then
    echo "--coverage-head-sha must identify the checked-out commit." >&2
    exit 1
  fi
  if [[ -n "${COVERAGE_BASE_SHA}" ]]; then
    if [[ ! "${COVERAGE_BASE_SHA}" =~ ^[0-9a-f]{40}$ ]] ||
       ! git -c "safe.directory=${GEOS_SRC_DIR}" -C "${GEOS_SRC_DIR}" \
             cat-file -e "${COVERAGE_BASE_SHA}^{commit}" 2>/dev/null; then
      echo "--coverage-base-sha must identify an available commit." >&2
      exit 1
    fi
  fi
  if [[ -n "${COVERAGE_BASELINE_SUMMARY}" ]]; then
    if [[ -z "${COVERAGE_BASE_SHA}" ||
          ! -f "${COVERAGE_BASELINE_SUMMARY}" ||
          -L "${COVERAGE_BASELINE_SUMMARY}" ]]; then
      echo "--coverage-baseline-summary requires a safe file and a base SHA." >&2
      exit 1
    fi
  fi

  # The pinned Clang CI image intentionally contains the compiler and TPLs but
  # not the compiler-native coverage reporter or profile runtime. Install the
  # matching LLVM-major packages in this ephemeral container.
  or_die apt-get update -qq
  or_die env DEBIAN_FRONTEND=noninteractive apt-get install -y \
               --no-install-recommends llvm-20 libclang-rt-20-dev
  export GEOS_CI_LLVM_PACKAGE_VERSIONS="$({
    dpkg-query -W -f='${Package}=${Version}\n' llvm-20 libclang-rt-20-dev
  } | LC_ALL=C sort)"
  or_die test -n "${GEOS_CI_LLVM_PACKAGE_VERSIONS}"

  for llvm_tool in llvm-cov-20 llvm-profdata-20 llvm-readelf-20; do
    or_die command -v "${llvm_tool}"
  done
  or_die mkdir -p "${COVERAGE_OUTPUT_DIR}"
  for coverage_artifact in \
    CMakeCache.txt \
    branch-summary.txt \
    coverage-objects.txt \
    coverage-status.md \
    coverage-summary.json \
    coverage-summary.md \
    ctest.log \
    index.html \
    llvm-summary.json \
    mapping-integrity.log \
    native-export.log \
    pr-coverage.json \
    pr-coverage.md \
    summary-export.log \
    toolchain.txt; do
    or_die rm -f -- "${COVERAGE_OUTPUT_DIR}/${coverage_artifact}"
  done
  or_die python3 -m unittest discover \
               -s "${GEOS_SRC_DIR}/scripts/tests" \
               -p 'test_*coverage*.py'
fi

# The -DBLT_MPI_COMMAND_APPEND="--allow-run-as-root;--oversubscribe" option is added for OpenMPI.
#
# OpenMPI prevents from running as `root` user by default.
# And by default user is `root` in a docker container.
# Using this option therefore offers a minimal and convenient way to run the unit tests.
#
# The option `--oversubscribe` tells OpenMPI to allow more MPI ranks than the node has cores.
# This is needed because our unit test `blt_mpi_smoke` is run in parallel with _hard coded_ 4 ranks.
# While some of our ci nodes may have less cores available.
#
# In case we have more powerful nodes, consider removing `--oversubscribe` and use `--use-hwthread-cpus` instead.
# This will tells OpenMPI to discover the number of hardware threads on the node,
# and use that as the number of slots available. (There is a distinction between threads and cores).
GEOS_BUILD_DIR=/tmp/geos-build
LLVM_COVERAGE_CMAKE_ARGS=(
  -DENABLE_COVERAGE=OFF
  -DGEOS_ENABLE_LLVM_SOURCE_COVERAGE=$([[ "${LLVM_SOURCE_COVERAGE}" = true ]] && echo ON || echo OFF)
)
if [[ "${LLVM_SOURCE_COVERAGE}" = true ]]; then
  LLVM_COVERAGE_CMAKE_ARGS+=(
    -DGEOS_BUILD_SHARED_LIBS=ON
    -DENABLE_DOXYGEN=OFF
  )
  LLVM_BUILD_PROFILE_DIR=/tmp/geos-coverage-build-profiles
  or_die mkdir -p "${LLVM_BUILD_PROFILE_DIR}"
  export LLVM_PROFILE_FILE="${LLVM_BUILD_PROFILE_DIR}/%p-%m.profraw"
fi

or_die python3 scripts/config-build.py \
               -hc ${HOST_CONFIG} \
               -bt ${CMAKE_BUILD_TYPE} \
               -bp ${GEOS_BUILD_DIR} \
               -ip ${GEOS_DIR} \
               ${BUILD_GENERATOR} \
               -DBLT_MPI_COMMAND_APPEND='"--allow-run-as-root;--oversubscribe"' \
               -DGEOS_INSTALL_SCHEMA=${GEOS_INSTALL_SCHEMA} \
               -DENABLE_HYPRE=${ENABLE_HYPRE} \
               -DENABLE_HYPREDRV=${ENABLE_HYPREDRV} \
               -DENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE} \
               -DENABLE_TRILINOS=${ENABLE_TRILINOS} \
               -DGEOS_LA_INTERFACE:PATH=${GEOS_LA_INTERFACE} \
               -DGEOS_ENABLE_BOUNDS_CHECK=${GEOS_ENABLE_BOUNDS_CHECK} \
               "${CMAKE_CUDA_ARCHITECTURES_ARGS[@]}" \
               "${CMAKE_NATIVE_ARCHITECTURE_ARGS[@]}" \
               ${SCCACHE_CMAKE_ARGS} \
               "${LLVM_COVERAGE_CMAKE_ARGS[@]}" \
               "${ATS_CMAKE_ARGS[@]}"

# The configuration step is now over, we can now move to the build directory for the build!
or_die cd ${GEOS_BUILD_DIR}

# Code style check
if [[ "${TEST_CODE_STYLE}" = true ]]; then
  or_die ctest --output-on-failure -R "testUncrustifyCheck"
  exit 0
fi

# Documentation check
if [[ "${TEST_DOCUMENTATION}" = true ]]; then
  or_die ctest --output-on-failure -R "testDoxygenCheck"
  exit 0
fi

# Performing the requested build.
if [[ "${BUILD_EXE_ONLY}" = true ]]; then
  or_die cmake --build . --parallel "${NPROC}" --target geosx
else
  or_die cmake --build . --parallel "${NPROC}"
  #or_die cmake --install .

  if [[ ! -z "${DATA_BASENAME_WE}" ]]; then
    # Here we pack the installation.
    # The `--transform` parameter provides consistency between the tarball name and the unpacked folder.
    echo "DATA_EXCHANGE_DIR=${DATA_EXCHANGE_DIR}"
    echo "DATA_BASENAME_WE=${DATA_BASENAME_WE}"
    echo "GEOS_TPL_DIR=${GEOS_TPL_DIR}"
    echo "GEOSX_TPL_DIR=${GEOSX_TPL_DIR}"
    GEOS_TPL_DIR=${GEOSX_TPL_DIR}
    echo tar czf ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar.gz --directory=${GEOS_TPL_DIR}/.. --transform "s|^./|${DATA_BASENAME_WE}/|" .
    or_die tar czf ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar.gz --directory=${GEOS_TPL_DIR}/.. --transform "s|^./|${DATA_BASENAME_WE}/|" .
  fi
fi

if [[ -n "${SCCACHE_BIN}" ]]; then
  echo "Capturing sccache post-build stats"
  SCCACHE_STATS_FILE="${GEOS_SRC_DIR}/.sccache-runtime/stats.txt"
  or_die mkdir -p "$(dirname "${SCCACHE_STATS_FILE}")"
  ${SCCACHE_BIN} --show-adv-stats > "${SCCACHE_STATS_FILE}"
  SCCACHE_STATS_STATUS=$?
  if [[ ${SCCACHE_STATS_STATUS} != 0 ]]; then
    echo ERROR ${SCCACHE_STATS_STATUS} command: ${SCCACHE_BIN} --show-adv-stats
    exit ${SCCACHE_STATS_STATUS}
  fi
fi

if [[ "${LLVM_SOURCE_COVERAGE}" = true ]]; then
  LLVM_TEST_PROFILE_DIR=/tmp/geos-coverage-profiles
  LLVM_REPORT_DIR=/tmp/geos-coverage-report
  LLVM_CTEST_PARALLEL="${CTEST_PARALLEL_LEVEL_ARG:-12}"
  or_die mkdir -p "${LLVM_TEST_PROFILE_DIR}" "${LLVM_REPORT_DIR}"

  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export LLVM_PROFILE_FILE="${LLVM_TEST_PROFILE_DIR}/%p-%m.profraw"

  # The CI container runs as root so it can install the compiler-matched LLVM
  # runtime. Open MPI requires this explicit, two-part acknowledgement before
  # it will launch the three MPI coverage smokes.
  if [[ "$(id -u)" -eq 0 ]]; then
    export OMPI_ALLOW_RUN_AS_ROOT=1
    export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
  fi

  {
    printf 'container image: %s\n' "${GEOS_CI_CONTAINER_IMAGE:-unknown}"
    printf 'container image ID: %s\n' "${GEOS_CI_CONTAINER_IMAGE_ID:-unknown}"
    printf 'container image digests: %s\n' "${GEOS_CI_CONTAINER_IMAGE_DIGESTS:-unknown}"
    or_die /usr/bin/clang-20 --version
    or_die llvm-cov-20 --version
    or_die llvm-profdata-20 --version
    or_die llvm-readelf-20 --version
    or_die cmake --version
    or_die dpkg-query -W -f='${Package}=${Version}\n' llvm-20 libclang-rt-20-dev
  } > "${COVERAGE_OUTPUT_DIR}/toolchain.txt"
  or_die test -s "${COVERAGE_OUTPUT_DIR}/toolchain.txt"

  coverage_ctest_status=0
  ctest --test-dir "${GEOS_BUILD_DIR}" \
        --parallel "${LLVM_CTEST_PARALLEL}" \
        --progress \
        --output-on-failure \
        --output-log "${COVERAGE_OUTPUT_DIR}/ctest.log" \
        -L '^coverage_smoke$' || coverage_ctest_status=$?

  coverage_report_status=0
  "${GEOS_SRC_DIR}/scripts/llvm_source_branch_coverage.sh" \
    "${GEOS_BUILD_DIR}" \
    "${LLVM_TEST_PROFILE_DIR}" \
    "${LLVM_REPORT_DIR}" || coverage_report_status=$?

  coverage_pr_status=0
  coverage_pr_ran=false
  if [[ ${coverage_report_status} -eq 0 && -n "${COVERAGE_BASE_SHA}" ]]; then
    coverage_pr_ran=true
    coverage_pr_args=(
      --repository "${GEOS_SRC_DIR}"
      --base-sha "${COVERAGE_BASE_SHA}"
      --head-sha "${COVERAGE_HEAD_SHA}"
      --native-info "${LLVM_REPORT_DIR}/native.info"
      --candidate-summary "${LLVM_REPORT_DIR}/coverage-summary.json"
      --output-json "${COVERAGE_OUTPUT_DIR}/pr-coverage.json"
      --output-markdown "${COVERAGE_OUTPUT_DIR}/pr-coverage.md"
    )
    if [[ -n "${COVERAGE_BASELINE_SUMMARY}" ]]; then
      coverage_pr_args+=(--baseline-summary "${COVERAGE_BASELINE_SUMMARY}")
    fi
    python3 "${GEOS_SRC_DIR}/scripts/compare_pr_coverage.py" \
      "${coverage_pr_args[@]}" || coverage_pr_status=$?
    if [[ ${coverage_pr_status} -eq 0 ]]; then
      echo "LLVM pull-request coverage report written to ${COVERAGE_OUTPUT_DIR}/pr-coverage.md"
    fi
  fi

  coverage_gate_status=0
  coverage_gate_ran=false
  if [[ ${coverage_report_status} -eq 0 ]]; then
    coverage_gate_ran=true
    python3 "${GEOS_SRC_DIR}/scripts/check_coverage_thresholds.py" \
      "${LLVM_REPORT_DIR}/coverage-summary.json" \
      "${GEOS_SRC_DIR}/.github/coverage-thresholds.json" \
      --markdown "${COVERAGE_OUTPUT_DIR}/coverage-summary.md" \
      || coverage_gate_status=$?
    if [[ ${coverage_gate_status} -eq 0 ]]; then
      echo "LLVM coverage policy report written to ${COVERAGE_OUTPUT_DIR}/coverage-summary.md"
    fi
    python3 "${GEOS_SRC_DIR}/scripts/render_coverage_review_summary.py" \
      --coverage-summary "${LLVM_REPORT_DIR}/coverage-summary.json" \
      --thresholds "${GEOS_SRC_DIR}/.github/coverage-thresholds.json" \
      --console || echo "LLVM coverage console summary unavailable." >&2
  fi

  for report_file in \
    branch-summary.txt \
    coverage-summary.json \
    coverage-objects.txt \
    llvm-summary.json \
    mapping-integrity.log \
    native-export.log \
    summary-export.log; do
    if [[ -f "${LLVM_REPORT_DIR}/${report_file}" ]]; then
      cp -- "${LLVM_REPORT_DIR}/${report_file}" "${COVERAGE_OUTPUT_DIR}/" \
        || coverage_report_status=1
    fi
  done
  cp -- "${GEOS_BUILD_DIR}/CMakeCache.txt" "${COVERAGE_OUTPUT_DIR}/" \
    || coverage_report_status=1

  coverage_overall_status=PASS
  if [[ ${coverage_ctest_status} -ne 0 ||
        ${coverage_report_status} -ne 0 ||
        ${coverage_gate_status} -ne 0 ||
        ${coverage_pr_status} -ne 0 ]]; then
    coverage_overall_status=FAIL
  fi

  coverage_phase_status()
  {
    if [[ "${2:-true}" != true ]]; then
      printf '⏭️ Not run'
      return
    fi
    if [[ "$1" -eq 0 ]]; then
      printf '✅ Pass'
    else
      printf '❌ Fail (%d)' "$1"
    fi
  }

  coverage_ctest_detail="See the retained CTest log."
  coverage_ctest_summary="$(
    grep -E '^[0-9]+% tests passed, [0-9]+ tests failed out of [0-9]+$' \
      "${COVERAGE_OUTPUT_DIR}/ctest.log" | tail -n 1 || true
  )"
  coverage_ctest_elapsed_line="$(
    grep -E '^Total Test time \(real\) = +[0-9.]+ sec$' \
      "${COVERAGE_OUTPUT_DIR}/ctest.log" | tail -n 1 || true
  )"
  coverage_ctest_summary_regex='^[0-9]+% tests passed, ([0-9]+) tests failed out of ([0-9]+)$'
  coverage_ctest_elapsed_regex='^Total Test time \(real\) = +([0-9.]+) sec$'
  if [[ "${coverage_ctest_summary}" =~ ${coverage_ctest_summary_regex} ]]; then
    coverage_ctest_failed="${BASH_REMATCH[1]}"
    coverage_ctest_total="${BASH_REMATCH[2]}"
    coverage_ctest_passed=$(( coverage_ctest_total - coverage_ctest_failed ))
    coverage_ctest_detail="${coverage_ctest_passed}/${coverage_ctest_total} passed"
    if (( coverage_ctest_failed > 0 )); then
      coverage_ctest_detail+=", ${coverage_ctest_failed} failed"
    fi
    if [[ "${coverage_ctest_elapsed_line}" =~ ${coverage_ctest_elapsed_regex} ]]; then
      coverage_ctest_detail+=" in ${BASH_REMATCH[1]} s"
    fi
    coverage_ctest_detail+=" with ${LLVM_CTEST_PARALLEL}-way scheduling"
  fi

  if [[ ${coverage_report_status} -eq 0 ]]; then
    coverage_report_detail="Profile mappings and report inputs validated"
  else
    coverage_report_detail="Report generation or artifact staging failed; inspect the logs"
  fi
  if [[ "${coverage_gate_ran}" != true ]]; then
    coverage_gate_detail="Skipped because the LLVM report was unavailable"
  elif [[ ${coverage_gate_status} -eq 0 ]]; then
    coverage_gate_detail="Exact repository thresholds evaluated below"
  else
    coverage_gate_detail="Repository threshold check failed"
  fi
  if [[ "${coverage_pr_ran}" != true ]]; then
    if [[ -n "${COVERAGE_BASE_SHA}" ]]; then
      coverage_pr_detail="Skipped because the LLVM report was unavailable"
    else
      coverage_pr_detail="Not a pull-request coverage run"
    fi
  elif [[ ${coverage_pr_status} -eq 0 ]]; then
    coverage_pr_detail="Changed-code coverage and exact-base comparison rendered below"
  else
    coverage_pr_detail="Patch-specific coverage analysis failed"
  fi

  coverage_status_markdown="${COVERAGE_OUTPUT_DIR}/coverage-status.md"
  coverage_status_tmp="${coverage_status_markdown}.tmp"
  coverage_image_id="${GEOS_CI_CONTAINER_IMAGE_ID:-unknown}"
  coverage_image_id_display="${coverage_image_id}"
  if (( ${#coverage_image_id_display} > 20 )); then
    coverage_image_id_display="${coverage_image_id_display:0:20}…"
  fi
  coverage_clang_version="$(
    /usr/bin/clang-20 --version | sed -nE \
      '1s/.*version[[:space:]]+([0-9]+([.][0-9]+)+).*/\1/p'
  )"
  if [[ -z "${coverage_clang_version}" ]]; then
    coverage_clang_version=20
  fi
  coverage_llvm_version="$(
    llvm-cov-20 --version | sed -nE \
      '1s/.*version[[:space:]]+([0-9]+([.][0-9]+)+).*/\1/p'
  )"
  if [[ -z "${coverage_llvm_version}" ]]; then
    coverage_llvm_version=20
  fi
  coverage_reviewer_args=(
    --coverage-summary "${LLVM_REPORT_DIR}/coverage-summary.json"
    --thresholds "${GEOS_SRC_DIR}/.github/coverage-thresholds.json"
    --tests-result "$(coverage_phase_status "${coverage_ctest_status}")"
    --tests-detail "${coverage_ctest_detail}"
  )
  if [[ "${coverage_pr_ran}" = true &&
        -f "${COVERAGE_OUTPUT_DIR}/pr-coverage.json" &&
        ! -L "${COVERAGE_OUTPUT_DIR}/pr-coverage.json" ]]; then
    coverage_reviewer_args+=(--pr-report "${COVERAGE_OUTPUT_DIR}/pr-coverage.json")
  fi
  coverage_reviewer_summary_status=0
  coverage_reviewer_summary="$(
    python3 "${GEOS_SRC_DIR}/scripts/render_coverage_review_summary.py" \
      "${coverage_reviewer_args[@]}"
  )" || coverage_reviewer_summary_status=$?
  if [[ ${coverage_reviewer_summary_status} -ne 0 ]]; then
    coverage_reviewer_summary=$'### Reviewer summary\n\nCoverage reviewer summary unavailable; inspect the detailed reports below.'
  fi
  {
    if [[ "${coverage_overall_status}" = PASS ]]; then
      printf '### ✅ LLVM source coverage passed\n\n'
    else
      printf '### ❌ LLVM source coverage failed\n\n'
    fi
    printf '%s\n' "${coverage_reviewer_summary}"
    printf '### Validation details\n\n'
    printf '| Validation | Result | Details |\n'
    printf '|---|:---:|---|\n'
    printf '| Coverage smoke suite | %s | %s |\n' \
      "$(coverage_phase_status "${coverage_ctest_status}")" "${coverage_ctest_detail}"
    printf '| LLVM report integrity | %s | %s |\n' \
      "$(coverage_phase_status "${coverage_report_status}")" "${coverage_report_detail}"
    printf '| Repository coverage policy | %s | %s |\n' \
      "$(coverage_phase_status "${coverage_gate_status}" "${coverage_gate_ran}")" \
      "${coverage_gate_detail}"
    printf '| PR coverage analysis | %s | %s |\n' \
      "$(coverage_phase_status "${coverage_pr_status}" "${coverage_pr_ran}")" \
      "${coverage_pr_detail}"
    printf '\n'
    printf '### Run configuration\n\n'
    printf '| Item | Value |\n'
    printf '|---|---|\n'
    printf '| Runner | `%s` |\n' "${GEOS_CI_RUNNER_LABEL:-unknown}"
    printf '| Container | `%s` |\n' "${GEOS_CI_CONTAINER_IMAGE:-unknown}"
    printf '| Immutable image ID (short) | `%s` |\n' "${coverage_image_id_display}"
    printf '| Toolchain | `Clang %s / LLVM %s` |\n' \
      "${coverage_clang_version}" "${coverage_llvm_version}"
    printf '| Build | `%s`, %s-way compilation |\n' "${CMAKE_BUILD_TYPE}" "${NPROC}"
    printf '| Test scheduling | %s CTest slots |\n' "${LLVM_CTEST_PARALLEL}"
    printf '\n'
  } > "${coverage_status_tmp}"
  or_die mv -- "${coverage_status_tmp}" "${coverage_status_markdown}"

  coverage_html_status=0
  python3 "${GEOS_SRC_DIR}/scripts/render_coverage_report_html.py" \
    --artifact-dir "${COVERAGE_OUTPUT_DIR}" \
    --output "${COVERAGE_OUTPUT_DIR}/index.html" \
    || coverage_html_status=$?
  if [[ ${coverage_html_status} -eq 0 ]]; then
    echo "LLVM source coverage HTML report written to ${COVERAGE_OUTPUT_DIR}/index.html"
  else
    echo "LLVM source coverage HTML report generation failed." >&2
  fi

  if [[ "${coverage_overall_status}" = FAIL ]]; then
    echo "LLVM source coverage failed: ctest=${coverage_ctest_status}, " \
         "report=${coverage_report_status}, gate=${coverage_gate_status}, " \
         "pr=${coverage_pr_status}, html=${coverage_html_status}." >&2
    exit 1
  fi
  if [[ ${coverage_html_status} -ne 0 ]]; then
    exit 1
  fi
  exit 0
fi

# Run the unit tests (excluding previously ran checks).
if [[ "${RUN_UNIT_TESTS}" = true ]]; then
  export OMP_NUM_THREADS=1
  echo "Running unit tests with OMP_NUM_THREADS=${OMP_NUM_THREADS}."
  if [ ${HOSTNAME} == 'streak.llnl.gov' ] || [ ${HOSTNAME} == 'streak2.llnl.gov' ]; then
    or_die ctest --output-on-failure --parallel -E "testUncrustifyCheck|testDoxygenCheck|testExternalSolvers"
  else
    or_die ctest --output-on-failure --parallel -E "testUncrustifyCheck|testDoxygenCheck"
  fi
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  PYTHON_MINOR_VERSION="$(${ATS_PYTHON_HOME}/bin/python3 -c 'import sys; print(sys.version_info.minor if sys.version_info.major == 3 else 99)')"
  if (( PYTHON_MINOR_VERSION < 12 )); then
    # fix the setuptools/distutils conflict
    export SETUPTOOLS_USE_DISTUTILS=stdlib
  else
    unset SETUPTOOLS_USE_DISTUTILS
  fi

  # We split the process in two steps. First installing the environment, then running the tests.
  or_die cmake --build . --target ats_environment
  # The system-site virtualenv inherits Ubuntu's NumPy 1.x-built Matplotlib,
  # while the GEOS Python packages require NumPy 2.x.
  or_die "${ATS_PYTHON_HOME}/bin/python3" -m pip install --disable-pip-version-check --upgrade matplotlib

  # The tests are not run using cmake (`cmake --build . --verbose  --target ats_run`)
  # because with ninja it swallows the output while all the
  # simulations are running.
  # We directly use the script instead...
  echo "Available baselines:"
  ls -lR /tmp/geos/baselines

  # This CI route configures geos-ats' Open MPI machine explicitly. Capture the
  # selected launcher's identity once, and fail before tests if an image change
  # would make the Open MPI placement environment below ineffective.
  echo "Diagnostic: mpirun --version"
  MPIRUN_VERSION_OUTPUT="$(/usr/bin/mpirun --version 2>&1)"
  MPIRUN_VERSION_STATUS=$?
  printf '%s\n' "${MPIRUN_VERSION_OUTPUT}"
  echo "Diagnostic 'mpirun --version' exit status: ${MPIRUN_VERSION_STATUS}"
  if [[ "${MPIRUN_VERSION_STATUS}" -ne 0 ]]; then
    echo "Open MPI integrated tests require /usr/bin/mpirun --version to succeed."
    exit "${MPIRUN_VERSION_STATUS}"
  fi
  if [[ "${MPIRUN_VERSION_OUTPUT}" != *"Open MPI"* ]]; then
    echo "Open MPI integrated tests require /usr/bin/mpirun to identify itself as Open MPI."
    exit 1
  fi

  # ATS accounts for the rank capacity of concurrent tests but does not assign
  # disjoint CPU affinities to their independent mpirun processes. Disable each
  # launcher's local rank binding so Linux can schedule ranks across the
  # container's allowed CPU set; this policy does not guarantee exclusive cores.
  export OMPI_MCA_hwloc_base_binding_policy=none
  echo "Open MPI rank binding disabled with OMPI_MCA_hwloc_base_binding_policy=${OMPI_MCA_hwloc_base_binding_policy}."

  # Report the resulting Open MPI placement policy. Each report is written to
  # the individual launch's stderr, which ATS retains in TestResults for the
  # packed diagnostic artifact.
  export OMPI_MCA_hwloc_base_report_bindings=1
  echo "Open MPI binding reports enabled with OMPI_MCA_hwloc_base_report_bindings=${OMPI_MCA_hwloc_base_report_bindings}."

  echo "Running integrated tests..."
  integratedTests/geos_ats.sh --baselineCacheDirectory /tmp/geos/baselines
  ATS_RUN_STATUS=$?

  PROCESS_LOGS_STATUS=0
  echo "Processing logs..."
  bin/geos_ats_process_tests_fails --directory integratedTests/TestResults &> integratedTests/TestResults/processedTestsLogs.txt || PROCESS_LOGS_STATUS=$?
  if [[ "${PROCESS_LOGS_STATUS}" -eq 0 ]]; then
    echo "Packing logs..."
    tar -czf ${DATA_EXCHANGE_DIR}/test_logs_${DATA_BASENAME_WE}.tar.gz integratedTests/TestResults || PROCESS_LOGS_STATUS=$?
  fi

  echo "Checking results..."
  bin/geos_ats_log_check integratedTests/TestResults/test_results.ini -y ${GEOS_SRC_DIR}/.integrated_tests.yaml &> $tempdir/log_check.txt
  LOG_CHECK_STATUS=$?
  cat $tempdir/log_check.txt

  if [[ "${ATS_RUN_STATUS}" -eq 0 && "${LOG_CHECK_STATUS}" -eq 0 ]] && grep -q "Overall status: PASSED" "$tempdir/log_check.txt"; then
    echo "IntegratedTests passed. No rebaseline required."
    INTEGRATED_TEST_EXIT_STATUS=0
  else
    echo "IntegratedTests failed. Rebaseline is required."

    # Rebaseline and pack into an archive
    echo "Rebaselining..."
    REBASELINE_STATUS=0
    integratedTests/geos_ats.sh -a rebaselinefailed || REBASELINE_STATUS=$?

    if [[ "${REBASELINE_STATUS}" -eq 0 ]]; then
      echo "Packing baselines..."
      integratedTests/geos_ats.sh -a pack_baselines --baselineArchiveName ${DATA_EXCHANGE_DIR}/baseline_${DATA_BASENAME_WE}.tar.gz --baselineCacheDirectory /tmp/geos/baselines || REBASELINE_STATUS=$?
    fi
    INTEGRATED_TEST_EXIT_STATUS=1
  fi

  echo "Done!"

  echo "The return code of the integrated tests is ${INTEGRATED_TEST_EXIT_STATUS}"
fi

# Cleaning the build directory.
or_die cmake --build . --target clean

# Clean the repository
or_die cd ${GEOS_SRC_DIR}/inputFiles
find . -name '*.pyc' -delete

# If we're here, either everything went OK or we have to deal with the integrated tests manually.
if [[ ! -z "${INTEGRATED_TEST_EXIT_STATUS+x}" ]]; then
  echo "Exiting the build process with exit status ${INTEGRATED_TEST_EXIT_STATUS} from the integrated tests."
  exit ${INTEGRATED_TEST_EXIT_STATUS}
else
  echo "Exiting the build process with exit status 0."
  exit 0
fi
