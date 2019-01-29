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
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    mkdir build-darwin-clang-debug;
    cd build-darwin-clang-debug;
else
    cd /home/geosx/geosx_repo
    export PATH=${PATHMOD}:$PATH
    or_die mkdir travis-build
    cd travis-build
fi

if [[ "$DO_BUILD" == "yes" ]] ; then
    or_die cmake \
           -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} \
           -DENABLE_MPI=ON -DMPI_C_COMPILER=${MPICC} -DMPI_CXX_COMPILER=${MPICXX} -DMPI_Fortran_COMPILER=${MPIFC} -DMPIEXEC=${MPIEXEC} -DMPIEXEC_EXECUTABLE=${MPIEXEC} \
           -DGEOSX_TPL_DIR=${GEOSX_TPL_DIR} \
           -DENABLE_SPHINX=OFF \
           ${CMAKE_EXTRA_FLAGS} ../src

#    if [[ ${CMAKE_EXTRA_FLAGS} == *COVERAGE* ]] ; then
#      or_die make -j 1
#    else
      or_die make -j 2 VERBOSE=1
#    fi
    if [[ "${DO_TEST}" == "yes" ]] ; then
      or_die ctest -V
    fi
fi

exit 0
