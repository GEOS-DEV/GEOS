#!/bin/bash
env

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

# Working in the root of the cloned repository
or_die cd $(dirname $0)/..

# The -DBLT_MPI_COMMAND_APPEND:STRING=--allow-run-as-root option is added for openmpi
# which prevents from running as root user by default.
# And by default, you are root in a docker container.
# Using this option therefore offers a minimal and convenient way
# to run the unit tests.
GEOSX_BUILD_DIR=/tmp/build
or_die python scripts/config-build.py \
              -hc host-configs/environment.cmake \
              -bt ${CMAKE_BUILD_TYPE} \
              -bp ${GEOSX_BUILD_DIR} \
              -ip ${GEOSX_DIR} \
              -DBLT_MPI_COMMAND_APPEND:STRING=--allow-run-as-root \
              -DENABLE_CUDA:BOOL=${ENABLE_CUDA:-OFF} \
              -DCMAKE_CUDA_FLAGS:STRING=\""${CMAKE_CUDA_FLAGS:-Unused}"\" \
              -DCUDA_TOOLKIT_ROOT_DIR:PATH=${CUDA_TOOLKIT_ROOT_DIR:-/usr/local/cuda} \
              -DCUDA_ARCH:STRING=${CUDA_ARCH:sm_70}

or_die cd ${GEOSX_BUILD_DIR}

# Code style check
if [[ "$*" == *--test-code-style* ]]; then
  or_die ctest -V -R "testUncrustifyCheck"
  exit 0
fi

# Documentation check
if [[ "$*" == *--test-documentation* ]]; then
  or_die ctest -V -R "testDoxygenCheck"
  exit 0
fi

or_die make -j $(nproc) VERBOSE=1
or_die make install VERBOSE=1

# Unit tests (excluding previously ran checks)
if [[ "$*" != *--disable-unit-tests* ]]; then
  or_die ctest -V -E "testUncrustifyCheck|testDoxygenCheck"
fi

exit 0
