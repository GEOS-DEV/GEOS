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
  --code-coverage
      run a code build and test.
  --ctest-parallel-level N
      Number of tests ctest may run in parallel.
  --data-basename output.tar.gz
      If some data needs to be extracted from the build, the argument will define the tarball. Has to be a `tar.gz`.
  --geos-enable-bounds-check
      Either ON or OFF (default is ON). Build geos with bounds check.
  --enable-hypre
      One of ON or OFF (default is ON). Build geos with hypre.
  --enable-hypredrv
      One of ON or OFF (default is OFF). Build geos with hypredrive support
      (requires hypre and a TPL image providing a hypredrive install).
      Overrides the TPL image host-config, which may set ENABLE_HYPREDRV ON.
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

args=$(or_die getopt -a -o h --long build-exe-only,cmake-build-type:,cmake-cuda-architectures:,code-coverage,ctest-parallel-level:,data-basename:,geos-enable-bounds-check:,enable-hypre:,enable-hypre-device:,enable-hypredrv:,enable-trilinos:,exchange-dir:,host-config:,install-dir-basename:,makefile,ninja,no-install-schema,no-run-unit-tests,nproc:,repository:,run-integrated-tests,sccache-config:,sccache-credentials:,test-code-style,test-documentation,use-native-architecture,use-sccache,help -- "$@")

# Variables with default values
BUILD_EXE_ONLY=false
BUILD_GENERATOR=""
GEOS_INSTALL_SCHEMA=true
HOST_CONFIG="host-configs/environment.cmake"
ENABLE_HYPRE=ON
ENABLE_HYPREDRV=OFF
ENABLE_HYPRE_DEVICE=CPU
GEOS_LA_INTERFACE=Hypre
RUN_UNIT_TESTS=true
RUN_INTEGRATED_TESTS=false
UPLOAD_TEST_BASELINES=false
TEST_CODE_STYLE=false
TEST_DOCUMENTATION=false
ENABLE_TRILINOS=OFF
CODE_COVERAGE=false
CTEST_PARALLEL_LEVEL_ARG=""
NPROC="$(nproc)"
GEOS_ENABLE_BOUNDS_CHECK=ON
SCCACHE_BIN=""
SCCACHE_CREDS=""
USE_SCCACHE=false
CMAKE_CUDA_ARCHITECTURES_ARGS=()
CMAKE_NATIVE_ARCHITECTURE_ARGS=()
ATS_CMAKE_ARGS=()
LCOV_CMAKE_ARGS=""

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
    --code-coverage)         CODE_COVERAGE=true;         shift;;
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

# Hypredrive support. GEOS_REQUIRE_HYPREDRV turns the silent
# "ENABLE_HYPREDRV forced OFF because HYPREDRV_DIR is missing" downgrade into a
# configuration error, so a CI job cannot silently test the legacy path instead.
HYPREDRV_CMAKE_ARGS=()
if [[ "${ENABLE_HYPREDRV}" = ON ]]; then
  HYPREDRV_CMAKE_ARGS=(-DHYPREDRV_DIR=${GEOSX_TPL_DIR}/hypredrive
                       -DGEOS_REQUIRE_HYPREDRV=ON)
fi

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
  ATS_ARGUMENTS="--machine openmpi --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_args=--allow-run-as-root --ats openmpi_procspernode=32 --ats openmpi_maxprocs=32 --ats cutoff=45m"
  ATS_CMAKE_ARGS=("-DATS_ARGUMENTS:STRING=\"${ATS_ARGUMENTS}\""
                  "-DPython3_ROOT_DIR=${ATS_PYTHON_HOME}"
                  "-DPython3_EXECUTABLE=${ATS_PYTHON_HOME}/bin/python3"
                  "-DATS_BASELINE_DIR=${ATS_BASELINE_DIR}"
                  "-DATS_WORKING_DIR=${ATS_WORKING_DIR}")
fi


if [[ "${CODE_COVERAGE}" = true ]]; then
  or_die apt-get update
  or_die apt-get install -y lcov

  LCOV_REAL=$(command -v lcov || true)
  if [[ -n "${LCOV_REAL}" ]]; then
    export GEOS_REAL_LCOV="${LCOV_REAL}"
    LCOV_WRAPPER=/tmp/geos-lcov-wrapper
    cat > "${LCOV_WRAPPER}" <<'EOF'
#!/bin/bash
set -e

extra_args=()
if "${GEOS_REAL_LCOV}" --version 2>&1 | grep -Eq 'LCOV version ([2-9]|[1-9][0-9])\.'; then
  for arg in "$@"; do
    case "${arg}" in
      --capture|-c)
        extra_args=(--ignore-errors mismatch,empty)
        break
        ;;
      --remove|-r)
        extra_args=(--ignore-errors unused)
        break
        ;;
    esac
  done
fi

exec "${GEOS_REAL_LCOV}" "${extra_args[@]}" "$@"
EOF
    or_die chmod +x "${LCOV_WRAPPER}"
    LCOV_CMAKE_ARGS="-DLCOV_EXECUTABLE=${LCOV_WRAPPER}"
  fi
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
               "${HYPREDRV_CMAKE_ARGS[@]}" \
               -DGEOS_LA_INTERFACE:PATH=${GEOS_LA_INTERFACE} \
               -DENABLE_COVERAGE=$([[ "${CODE_COVERAGE}" = true ]] && echo 1 || echo 0) \
               -DGEOS_ENABLE_BOUNDS_CHECK=${GEOS_ENABLE_BOUNDS_CHECK} \
               "${CMAKE_CUDA_ARCHITECTURES_ARGS[@]}" \
               "${CMAKE_NATIVE_ARCHITECTURE_ARGS[@]}" \
               ${SCCACHE_CMAKE_ARGS} \
               ${LCOV_CMAKE_ARGS} \
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
  or_die cmake --build . -j $NPROC --target geosx
else
  or_die cmake --build . -j $NPROC
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

if [[ "${CODE_COVERAGE}" = true ]]; then
  export OMP_NUM_THREADS=1
  or_die cmake --build . --target coreComponents_coverage
  or_die cp -r ${GEOS_BUILD_DIR}/coreComponents_coverage.info.cleaned ${GEOS_SRC_DIR}/geos_coverage.info.cleaned
fi

# Run the unit tests (excluding previously ran checks).
if [[ "${RUN_UNIT_TESTS}" = true ]]; then
  if [ ${HOSTNAME} == 'streak.llnl.gov' ] || [ ${HOSTNAME} == 'streak2.llnl.gov' ]; then
    or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck|testExternalSolvers"
  else
    or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck"
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

  echo "Running integrated tests..."
  integratedTests/geos_ats.sh --baselineCacheDirectory /tmp/geos/baselines
  ATS_RUN_STATUS=$?

  if [[ "${ENABLE_HYPREDRV}" = ON ]]; then
    # On a hypredrive-enabled build the previous pass exercised the hypredrive
    # solver path. Verify it actually ran (guards against silent fallback), then
    # re-run the suite through the legacy hypre path against the same baselines,
    # so hypredrive-vs-legacy equivalence is checked within a single job.
    HYPREDRV_BANNER_STATUS=0
    grep -rl "hypredrive input" integratedTests/workingDir > /dev/null 2>&1 || HYPREDRV_BANNER_STATUS=$?
    if [[ "${HYPREDRV_BANNER_STATUS}" -ne 0 ]]; then
      echo "ERROR: no 'hypredrive input' banner found in any integrated-test log; the hypredrive path was not exercised."
    fi

    echo "Re-running integrated tests through the legacy hypre path..."
    # Keep the hypredrive pass results; the second pass overwrites TestResults.
    cp integratedTests/TestResults/test_results.ini $tempdir/test_results_hypredrive.ini
    python3 ${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py harvest \
            integratedTests/workingDir -o $tempdir/iterations_hypredrive.json
    integratedTests/geos_ats.sh -a veryclean
    GEOS_HYPREDRV_FORCE_LEGACY=1 integratedTests/geos_ats.sh --baselineCacheDirectory /tmp/geos/baselines
    ATS_LEGACY_RUN_STATUS=$?

    # Restart checks compare solution fields only; additionally require iteration-count
    # parity between the two passes so solver-quality regressions cannot pass silently.
    python3 ${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py harvest \
            integratedTests/workingDir -o $tempdir/iterations_legacy.json
    ITER_PARITY_STATUS=0
    python3 ${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py compare \
            $tempdir/iterations_hypredrive.json $tempdir/iterations_legacy.json || ITER_PARITY_STATUS=$?
  fi

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
  if [[ -f $tempdir/test_results_hypredrive.ini ]]; then
    # Also gate on the hypredrive pass saved before the forced-legacy re-run.
    bin/geos_ats_log_check $tempdir/test_results_hypredrive.ini -y ${GEOS_SRC_DIR}/.integrated_tests.yaml &> $tempdir/log_check_hypredrive.txt
    cat $tempdir/log_check_hypredrive.txt
    if ! grep -q "Overall status: PASSED" "$tempdir/log_check_hypredrive.txt"; then
      echo "Hypredrive-path integrated tests did not pass."
      LOG_CHECK_STATUS=1
      # Force the rebaseline branch below to report failure.
      echo "Overall status: FAILED (hypredrive pass)" >> $tempdir/log_check.txt
      sed -i 's/Overall status: PASSED//' $tempdir/log_check.txt
    fi
  fi

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
