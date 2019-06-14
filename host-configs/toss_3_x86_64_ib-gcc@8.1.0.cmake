set(CONFIG_NAME "toss_3_x86_64_ib-gcc@8.1.0" CACHE PATH "") 

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gfortran" CACHE PATH "")

# set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "")
# set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "")
# set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set(GEOSX_TPL_ROOT_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/" CACHE PATH "")
set(GEOSX_TPL_DIR "${GEOSX_TPL_ROOT_DIR}/install-${CONFIG_NAME}-release" CACHE PATH "")

set(SPHINX_EXECUTABLE "/collab/usr/gapps/python/build/spack-toss3.2/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.14-7rci3jkmuht2uiwp433afigveuf4ocnu/bin/sphinx-build" CACHE PATH "" FORCE)
set(UNCRUSTIFY_EXECUTABLE "${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify" CACHE PATH "" FORCE)
set(DOXYGEN_EXECUTABLE "/usr/gapps/GEOS/geosx/thirdPartyLibs/doxygen/bin/doxygen" CACHE PATH "" FORCE)

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

option(RAJA_ENABLE_TBB "" OFF)

option(ENABLE_CALIPER "Enables CALIPER" ON)

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "ON"        CACHE BOOL "" FORCE)
set(RAJA_ENABLE_OPENMP "ON"        CACHE BOOL "" FORCE)

# Use mkl instead of the lapack suite
set(ENABLE_LAPACK_SUITE OFF CACHE BOOL "" FORCE)
set(ENABLE_MKL ON CACHE BOOL "" FORCE)
set(MKL_ROOT "/usr/tce/packages/mkl/mkl-2019.0" CACHE STRING "" FORCE)
set(MKL_LIBRARY_NAMES libmkl_rt.so CACHE STRING "" FORCE)
