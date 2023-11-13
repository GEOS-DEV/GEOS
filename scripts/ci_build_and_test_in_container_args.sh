#!/bin/bash
set -o pipefail

printenv

SCRIPT_NAME=$0
echo "Running CLI ${SCRIPT_NAME} $@"

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
  [ --cmake-build-type ... ]
  [ --data-basename ... ]
  [ --exchange-dir ... ]
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

# First working in the root of the cloned repository.
# Then we'll move to the build dir.
or_die cd $(dirname $0)/..

args=$(or_die getopt -a -o h --long build-exe-only,cmake-build-type:,data-basename:,exchange-dir:,no-run-unit-tests,host-config:,install-dir-basename:,no-install-schema,install-dir:,test-code-style,test-documentation,sccache-credentials:,repository:,run-integrated-tests,help -- "$@")
if [[ $? -gt 0 ]]; then
  echo "Error after getopt"
  echo "which getop"
  echo $(which getopt)
  usage
fi

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
    --build-exe-only)
      BUILD_EXE_ONLY=true
      RUN_UNIT_TESTS=false
      shift;;
    --cmake-build-type)      CMAKE_BUILD_TYPE=$2;        shift 2;;
    --data-basename)
      DATA_BASENAME=$2
      DATA_BASENAME_WE=${DATA_BASENAME%%.*}
      DATA_BASENAME_EXT=${DATA_BASENAME#*.}
      if [[ ${DATA_BASENAME_EXT} != "tar.gz" ]] ; then
          echo "The script ${SCRIPT_NAME} can only pack data into a '.tar.gz' file."
          exit 1
      fi
      unset DATA_BASENAME DATA_BASENAME_EXT  # We do not need those anymore
      shift 2;;
    --exchange-dir)          DATA_EXCHANGE_DIR=$2;       shift 2;;
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

if [[ -z "${GEOSX_DIR}" ]]; then
  echo "Installation folder undefined. Set to default value '/dev/null'. Please define it using '--install-dir-basename'."
  GEOSX_DIR=/dev/null
fi

if [[ ! -z "${SCCACHE_CREDS}" ]]; then
  or_die mkdir -p ${HOME}/.config/sccache
  or_die cat <<EOT >> ${HOME}/.config/sccache/config
[cache.gcs]
rw_mode = "READ_WRITE"
cred_path = "${GEOS_SRC_DIR}/${SCCACHE_CREDS}"
bucket = "geos-dev"
key_prefix = "sccache"
EOT

  or_die cat ${HOME}/.config/sccache/config

  SCCACHE_CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=${SCCACHE} -DCMAKE_CUDA_COMPILER_LAUNCHER=${SCCACHE}"

  echo "sccache initial state"
  ${SCCACHE} --show-stats
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  echo "We should be running the integrated tests."
  or_die apt-get install -y virtualenv python3-dev python-is-python3
  ATS_PYTHON_HOME=/tmp/run_integrated_tests_virtualenv
  or_die virtualenv ${ATS_PYTHON_HOME}
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

# The configuration step is done, we now move to the build dir for the build!
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
    # Here we pack the installation. 
    # The `--transform` parameter provides consistency between the tarball name and the unpacked folder.
    or_die tar czf ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar.gz --directory=${GEOSX_TPL_DIR}/.. --transform "s/^./${DATA_BASENAME_WE}/" .
  fi
fi

# Unit tests (excluding previously ran checks)
if [[ "${RUN_UNIT_TESTS}" = true ]]; then
  or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck"
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  # We split the process in two steps. First installing the environment, then runnint the tests.
  # The tests are not run using ninja because it swallows the output shile all the simulations are running.
  # We directly use the script instead.
  or_die ninja ats_environment
  # ninja --verbose ats_run
  cat /tmp/build/integratedTests/geos_ats.sh
  # integratedTests/geos_ats.sh --failIfTestsFail
  # /tmp/build/bin/run_geos_ats /tmp/build/bin --workingDir /tmp/geos/integratedTests/tests/allTests/simplePDE --logs /tmp/build/integratedTests/TestResults --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_args=--allow-run-as-root --ats openmpi_procspernode=2 --ats openmpi_maxprocs=2 --machine openmpi
  /tmp/build/bin/run_geos_ats /tmp/build/bin --workingDir /tmp/geos/integratedTests/tests/allTests/simplePDE --logs /tmp/build/integratedTests/TestResults --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_args=--allow-run-as-root --ats openmpi_procspernode=2 --ats openmpi_maxprocs=2 --machine openmpi --failIfTestsFail
  INTEGRATED_TEST_EXIT_STATUS=$?
  echo "The return code of the integrated tests is ${INTEGRATED_TEST_EXIT_STATUS}"

  # Whatever the result of the integrated tests, we want to pack both the logs and the computed results.
  # They are not in the same folder, so we do it in 2 steps.
  # The `--transform` parameter is here to separate the two informations (originally in a folder with the same name)
  # in two different folder with meaningful names when unpacking. 
  or_die tar cfM ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar --directory ${GEOS_SRC_DIR}    --transform "s/^integratedTests/${DATA_BASENAME_WE}\/repo/" integratedTests
  or_die tar rfM ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar --directory ${GEOSX_BUILD_DIR} --transform "s/^integratedTests/${DATA_BASENAME_WE}\/logs/" integratedTests
  or_die gzip ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar
fi

if [[ ! -z "${SCCACHE_CREDS}" ]]; then
  echo "sccache final state"
  or_die ${SCCACHE} --show-stats
fi

if [[ ! -z "${INTEGRATED_TEST_EXIT_STATUS+x}" ]]; then
  echo "Exiting the build process with exit status ${INTEGRATED_TEST_EXIT_STATUS} from the integrated tests."
  exit ${INTEGRATED_TEST_EXIT_STATUS}
else
  echo "Exiting the build process with exit status 0."
  exit 0
fi
