#!/bin/bash
#env
function or_die () {
    "$@"
    local status=$?
    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}
echo ${1:-0}

if [ "${1:-0}" = "1" ] ; then
    echo "local execution"
    export DO_BUILD="yes"
    export CMAKE_EXTRA_FLAGS=""
    export CC=/opt/local/bin/gcc-mp-7
    export CXX=/opt/local/bin/g++-mp-7
    export MPICC=/opt/local/bin/mpicc-mpich-gcc7
    export MPICXX=/opt/local/bin/mpicxx-mpich-gcc7
    rm -rf travis-build
else
    #    source ~/.bashrc
    cd /home/geosx/geosx_repo
#    [[ -d /opt/intel ]] && . /opt/intel/bin/compilervars.sh intel64
fi


pwd
echo $(ls)
or_die mkdir travis-build
cd travis-build
pwd
echo $(ls)
if [[ "$DO_BUILD" == "yes" ]] ; then
    or_die cmake \
           -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} \
           -DENABLE_MPI=ON -DMPI_C_COMPILER=${MPICC} -DMPI_CXX_COMPILER=${MPICXX} -DMPI_EXEC=mpirun \
           -DGEOSX_TPL_DIR=/home/geosx/thirdPartyLibs/install-default-release \
           -DSPHINX_EXECUTABLE=/usr/bin/sphinx-build \
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
