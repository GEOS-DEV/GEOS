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

cd /home/geosx/GEOSX
python scripts/config-build.py -hc hostconfigs/environment.cmake -bt ${CMAKE_BUILD_TYPE}




exit 0
