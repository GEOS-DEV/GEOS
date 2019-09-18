#!/bin/bash
env
function or_die () {
    "$@"
    local status=$?
    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

#export CMAKE_BUILD_TYPE=Debug
or_die cd /home/geosx/GEOSX
or_die python scripts/config-build.py -hc host-configs/environment.cmake -bt ${CMAKE_BUILD_TYPE} -DENABLE_GEOSX_PTP:BOOL=ON
or_die cd build-environment-*
or_die make -j 2 VERBOSE=1
or_die ctest -V


exit 0
