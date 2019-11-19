set(CONFIG_NAME "toss_3_x86_64_ib-clang@9.0.0" CACHE PATH "")

set(CMAKE_C_COMPILER /opt/moose/llvm-8.0.0/bin/clang CACHE PATH "")
set(CMAKE_CXX_COMPILER /opt/moose/llvm-8.0.0/bin/clang++ CACHE PATH "")
set(CMAKE_Fortran_COMPILER /opt/moose/gcc-9.2.0/bin/gfortran CACHE PATH "")

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "")
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             /opt/moose/mpich-3.3/gcc-9.2.0/ CACHE PATH "")
set(MPI_C_COMPILER       ${MPI_HOME}/bin/mpicc   CACHE PATH "")
set(MPI_CXX_COMPILER     ${MPI_HOME}/bin/mpicxx  CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

set(MPIEXEC              /usr/bin/srun CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(GEOSX_TPL_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs CACHE PATH "")
set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/11-11-19/install-${CONFIG_NAME}-release CACHE PATH "")

#set(SPHINX_EXECUTABLE /Home/projects/cmake-3.15.2/Utilities/Sphinx CACHE PATH "")
set(DOXYGEN_EXECUTABLE /usr/bin/doxygen CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)


set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI ON CACHE BOOL "")
#I don't know what papi is and I also don't know where it is so I will just comment this out
#set(PAPI_PREFIX /usr/tce/packages/papi/papi-5.4.3 CACHE PATH "")

set(ENABLE_OPENMP ON CACHE BOOL "")
set(CUDA_ENABLED OFF CACHE BOOL "")

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /opt/moose/miniconda/)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/libmkl_core.so
                  CACHE STRING "")


set( ENABLE_TOTALVIEW_OUTPUT ON CACHE BOOL "Enables Totalview custom view" FORCE)

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "")
