#!/bin/bash
# set -euo pipefail
set -eo pipefail

printenv

echo "Running CLI $0 $@"

echo "running nproc"
nproc

echo "running free -m"
free -m

echo "ls /tmp/geos"
ls /tmp/geos

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

function usage () {
>&2 cat << EOF
Usage: $0
  [ --build-exe-only ]
  [ --cmake-build-type ... ]
  [ --data-basename-we ... ]
  [ --exchange ... ]
  [ --host-config ... ]
  [ --no-install-schema ]
  [ --no-run-unit-tests ]
  [ --run-integrated-tests ]
  [ --sccache-credentials ... ]
  [ --test-code-style ]
  [ --test-documentation ]
  [ -h | --help ]
EOF
exit 1
}

# Working in the root of the cloned repository
or_die cd $(dirname $0)/..

args=$(getopt -a -o h --long build-exe-only,cmake-build-type:,data-basename-we:,exchange:,no-run-unit-tests,host-config:,install-dir-basename:,no-install-schema,install-dir:,test-code-style,test-documentation,sccache-credentials:,repository:,run-integrated-tests,help -- "$@")
if [[ $? -gt 0 ]]; then
  echo "Error after getopt"
  echo "which getop"
  echo $(which getopt)
  usage
fi
# TODO or_die args=$(...)

# Variables with default values
BUILD_EXE_ONLY=false
GEOSX_INSTALL_SCHEMA=true
HOST_CONFIG="host-configs/environment.cmake"
RUN_UNIT_TESTS=true
RUN_INTEGRATED_TESTS=false
TEST_CODE_STYLE=false
TEST_DOCUMENTATION=false

eval set -- ${args}
while :
do
  case $1 in
    --build-exe-only)        BUILD_EXE_ONLY=true; RUN_UNIT_TESTS=false; shift;;
    --cmake-build-type)      CMAKE_BUILD_TYPE=$2;        shift 2;;
    --data-basename-we)      DATA_BASENAME_WE=$2;        shift 2;;
    --exchange)              DATA_EXCHANGE_DIR=$2;       shift 2;;
    --host-config)           HOST_CONFIG=$2;             shift 2;;
    --install-dir-basename ) GEOSX_DIR=${GEOSX_TPL_DIR}/../$2; shift 2;;
    --no-install-schema)     GEOSX_INSTALL_SCHEMA=false; shift;;
    --no-run-unit-tests)     RUN_UNIT_TESTS=false;       shift;;
    --repository)            GEOS_SRC_DIR=$2;            shift 2;;
    --run-integrated-tests)  RUN_INTEGRATED_TESTS=true;  shift;;
    --sccache-credentials)   SCCACHE_CREDS=$2;           shift 2;;
    --test-code-style)       TEST_CODE_STYLE=true;       shift;;
    --test-documentation)    TEST_DOCUMENTATION=true;    shift;;
    -h | --help)             usage;                      shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

# if [[ $# -eq 0 ]]; then
#   usage
# fi

if [[ -z "${GEOS_SRC_DIR}" ]]; then
  echo "Variable GEOS_SRC_DIR is either empty or not defined. Please define it using '--repository'."
  exit 1
fi

# if [[ -z "${GEOSX_DIR}" ]]; then
#   echo "Variable GEOSX_DIR is either empty or not defined. Please define it using '--install-dir-basename'."
#   exit 1
# fi
# GEOSX_DIR=${GEOSX_TPL_DIR}/../GEOSX-${GITHUB_SHA:0:7}

if [[ ! -z "${SCCACHE_CREDS}" ]]; then
  mkdir -p ${HOME}/.config/sccache
  cat <<EOT >> ${HOME}/.config/sccache/config
[cache.gcs]
rw_mode = "READ_WRITE"
cred_path = "${GEOS_SRC_DIR}/${SCCACHE_CREDS}"
bucket = "geos-dev"
key_prefix = "sccache"
EOT

  cat ${HOME}/.config/sccache/config

  SCCACHE_CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=${SCCACHE} -DCMAKE_CUDA_COMPILER_LAUNCHER=${SCCACHE}"

  echo "sccache initial state"
  ${SCCACHE} --show-stats
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  echo "We should be running the integrated tests."
  apt-get install -y virtualenv python3-dev python-is-python3
  ATS_PYTHON_HOME=/tmp/run_integrated_tests_virtualenv
  virtualenv ${ATS_PYTHON_HOME}
  ATS_CMAKE_ARGS="-DATS_ARGUMENTS=\"--machine openmpi --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_args=--allow-run-as-root --ats openmpi_procspernode=2 --ats openmpi_maxprocs=2\" -DPython3_ROOT_DIR=${ATS_PYTHON_HOME}"
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
GEOSX_BUILD_DIR=/tmp/build
or_die python3 scripts/config-build.py \
               -hc ${HOST_CONFIG} \
               -bt ${CMAKE_BUILD_TYPE} \
               -bp ${GEOSX_BUILD_DIR} \
               -ip ${GEOSX_DIR} \
               --ninja \
               -DBLT_MPI_COMMAND_APPEND='"--allow-run-as-root;--oversubscribe"' \
               -DGEOSX_INSTALL_SCHEMA=${GEOSX_INSTALL_SCHEMA} \
               ${SCCACHE_CMAKE_ARGS} \
               ${ATS_CMAKE_ARGS}

or_die cd ${GEOSX_BUILD_DIR}

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

# "Make" target check (builds geosx executable target only if true)
# Use one process to prevent out-of-memory error
if [[ "${BUILD_EXE_ONLY}" = true ]]; then
  or_die ninja -j $(nproc) geosx
else
  or_die ninja -j $(nproc)
  or_die ninja install

  if [[ ! -z "${DATA_BASENAME_WE}" ]]; then
    # Here we pack the installation
   tar czf ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar.gz --directory=${GEOSX_TPL_DIR}/.. --transform "s/^/${DATA_BASENAME_WE}\//" .
  fi
fi

# Unit tests (excluding previously ran checks)
if [[ "${RUN_UNIT_TESTS}" = true ]]; then
  or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck"
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  echo "cwd is ${PWD}"
  # We split the process in two steps. First installing the environment, then runnint the tests.
  # The tests are not run using ninja because it swallows the output shile all the simulations are running.
  # We directly use the script instead.
  or_die ninja ats_environment
  # ninja --verbose ats_run
  cat /tmp/build/integratedTests/geos_ats.sh
  # integratedTests/geos_ats.sh --failIfTestsFail
  # integratedTests/geos_ats.sh
  /tmp/build/bin/run_geos_ats /tmp/build/bin --workingDir /tmp/geos/integratedTests/tests/allTests/simplePDE --logs /tmp/build/integratedTests/TestResults --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_args=--allow-run-as-root --ats openmpi_procspernode=2 --ats openmpi_maxprocs=2 --machine openmpi
  exit_status=$?
  echo "The return code is ${exit_status}"

  tar cfM ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar --directory ${GEOS_SRC_DIR}    --transform 's/^integratedTests/integratedTests\/repo/' integratedTests
  tar rfM ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar --directory ${GEOSX_BUILD_DIR} --transform 's/^integratedTests/integratedTests\/logs/' integratedTests
  gzip ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar
fi

if [[ ! -z "${SCCACHE_CREDS}" ]]; then
  echo "sccache final state"
  ${SCCACHE} --show-stats
fi
exit 0
