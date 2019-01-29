

set(CONFIG_NAME "toss_3_x86_64_ib-clang@6.0.0" CACHE PATH "")


set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-6.0.0/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-6.0.0/bin/clang++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/tce/packages/intel/intel-18.0.1/bin/ifort" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-6.0.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set( GEOSX_TPL_ROOT_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/" CACHE PATH "" )
set( GEOSX_TPL_DIR "${GEOSX_TPL_ROOT_DIR}/install-${CONFIG_NAME}-release" CACHE PATH "" )

set(SPHINX_EXECUTABLE "/collab/usr/gapps/python/build/spack-toss3.2/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.14-7rci3jkmuht2uiwp433afigveuf4ocnu/bin/sphinx-build" CACHE PATH "" FORCE)
set(UNCRUSTIFY_EXECUTABLE "${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify" CACHE PATH "" FORCE )
set(DOXYGEN_EXECUTABLE "${GEOSX_TPL_ROOT_DIR}/doxygen/bin/doxygen" CACHE PATH "" FORCE )

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
#set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )

option( RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" ON )
set(ENABLE_PAPI "ON" CACHE PATH "" FORCE)
set(PAPI_PREFIX "/usr/tce/packages/papi/papi-5.4.3" CACHE PATH "" FORCE)

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)

