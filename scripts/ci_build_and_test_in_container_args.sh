#!/bin/bash
set -euo pipefail

env

echo "running nproc"
nproc

echo "running free -m"
free -m

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
  [ --cmake-build-type ]
  [ --no-run-unit-tests ]
  [ --gcp-credential-file ]
  [ --host-config ]
  [ --no-install-schema ]
  [ --install-dir ]
  [ --test-code-style ]
  [ --test-documentation ]
  [ --no-use-sccache ]
  [ -h | --help ]
EOF
exit 1
}

# Working in the root of the cloned repository
or_die cd $(dirname $0)/..

# args=$(getopt -a -o h --long build-exe-only,cmake-build-type:,no-run-unit-tests,gcp-credential-file:,host-config:,no-install-schema,install-dir:,test-code-style,test-documentation,no-use-sccache,help -- "$@")
# if [[ $? -gt 0 ]]; then
#   usage
# fi
or_die args=$(getopt -a -o h --long build-exe-only,cmake-build-type:,no-run-unit-tests,gcp-credential-file:,host-config:,no-install-schema,install-dir:,test-code-style,test-documentation,no-use-sccache,help -- "$@")

# Variables and default values
BUILD_EXE_ONLY=0
CMAKE_BUILD_TYPE=unset
RUN_UNIT_TESTS=1
GCP_CREDENTIAL_FILE=unset
GEOSX_DIR=unset
GEOSX_INSTALL_SCHEMA=1
HOST_CONFIG=unset
TEST_CODE_STYLE=0
TEST_DOCUMENTATION=0
USE_SCCACHE=1

eval set -- ${args}
while :
do
  case $1 in
    --build-exe-only)      BUILD_EXE_ONLY=1;       shift;;
    --cmake-build-type)    CMAKE_BUILD_TYPE=$2;    shift 2;;
    --gcp-credential-file) GCP_CREDENTIAL_FILE=$2; shift 2;;
    --host-config)         HOST_CONFIG=$2;         shift 2;;
    --install-dir)         GEOSX_DIR=$2;           shift 2;;
    --no-install-schema)   GEOSX_INSTALL_SCHEMA=0; shift;;
    --no-run-unit-tests)   RUN_UNIT_TESTS=0;       shift;;
    --no-use-sccache)      USE_SCCACHE=0;          shift;;
    --test-code-style)     TEST_CODE_STYLE=1;      shift;;
    --test-documentation)  TEST_DOCUMENTATION=1;   shift;;
    -h | --help)           usage;                  shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

if [[ $# -eq 0 ]]; then
  usage
fi

if [[ -z "${HOST_CONFIG}" ]]; then
  echo "Variable \"HOST_CONFIG\" is undefined or empty. Define it using '--host-config'."
  exit 1
fi

if [[ -z "${CMAKE_BUILD_TYPE}" ]]; then
  echo "Variable \"CMAKE_BUILD_TYPE\" is undefined or empty. Define it using '--cmake-build-type'."
  exit 1
fi

if [[ -z "${GEOSX_DIR}" ]]; then
  echo "Variable \"GEOSX_DIR\" is undefined or empty. Define it using '--install-dir'."
  exit 1
fi

if [[ ${USE_SCCACHE} ]]; then
  mkdir -p ${HOME}/.config/sccache
  cat <<EOT >> ${HOME}/.config/sccache/config
[cache.gcs]
rw_mode = "READ_WRITE"
cred_path = "/opt/gcs/credentials.json"
bucket = "geos-dev"
key_prefix = "sccache"
EOT

  SCCACHE_CMAKE_PARAMETERS="-DCMAKE_CXX_COMPILER_LAUNCHER=${SCCACHE} -DCMAKE_CUDA_COMPILER_LAUNCHER=${SCCACHE}"

  echo "sccache initial state"
  ${SCCACHE} --show-stats
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
               -DGEOSX_INSTALL_SCHEMA=${GEOSX_INSTALL_SCHEMA} ${SCCACHE_CMAKE_PARAMETERS}

or_die cd ${GEOSX_BUILD_DIR}

# Code style check
# if [[ "$*" == *--test-code-style* ]]; then
if [[ ${TEST_CODE_STYLE} ]]; then
  or_die ctest --output-on-failure -R "testUncrustifyCheck"
  exit 0
fi

# Documentation check
if [[ ${TEST_DOCUMENTATION} ]]; then
  or_die ctest --output-on-failure -R "testDoxygenCheck"
  exit 0
fi

# "Make" target check (builds geosx executable target only if true)
# Use one process to prevent out-of-memory error
if [[ ${BUILD_EXE_ONLY} ]]; then
  or_die ninja -j $(nproc) geosx
else
  or_die ninja -j $(nproc)
  or_die ninja install
fi

# Unit tests (excluding previously ran checks)
if [[ ${RUN_UNIT_TESTS} ]]; then
  or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck"
fi
 
if [[ ${USE_SCCACHE} ]]; then
  echo "sccache final state"
  ${SCCACHE} --show-stats
fi

exit 0
