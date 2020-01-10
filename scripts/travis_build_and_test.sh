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

or_die cd /tmp/GEOSX
# The -DBLT_MPI_COMMAND_APPEND:STRING=--allow-run-as-root option is added for openmpi
# which prevents from running as root user by default.
# And by default, you are root in a docker container.
# Using this option therefore offers a minimal and convenient way
# to run the unit tests.
or_die python scripts/config-build.py \
              -hc host-configs/environment.cmake \
              -bt ${CMAKE_BUILD_TYPE} \
              -bp /tmp/build \
              -DGEOSX_TPL_DIR=$GEOSX_TPL_DIR \
              -DENABLE_GEOSX_PTP:BOOL=ON \
              -DBLT_MPI_COMMAND_APPEND:STRING=--allow-run-as-root \
              -DENABLE_CUDA:BOOL=${ENABLE_CUDA:-OFF} \
              -DCMAKE_CUDA_FLAGS:STRING=\""${CMAKE_CUDA_FLAGS:-Unused}"\"
              
or_die cd /tmp/build
or_die make -j $(nproc) VERBOSE=1

# Unit tests
if [[ "$*" != *-disable-unit-tests* ]]; then
  or_die ctest -V
fi 

exit 0
