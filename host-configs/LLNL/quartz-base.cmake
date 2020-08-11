###############################################################################
#
# Base configuration for LC builds
# Calling configuration file must define the following CMAKE variables:
#
# CONFIG_NAME
# CMAKE_C_COMPILER
# CMAKE_CXX_COMPILER
# CMAKE_Fortran_COMPILER
# MPI_HOME
# 
#
###############################################################################

set(ENABLE_FORTRAN OFF CACHE BOOL "")
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER       ${MPI_HOME}/bin/mpicc   CACHE PATH "")
set(MPI_CXX_COMPILER     ${MPI_HOME}/bin/mpicxx  CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

set(MPIEXEC              /usr/bin/srun CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(GEOSX_TPL_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs CACHE PATH "")
set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/2020-07-08/install-${CONFIG_NAME}-release CACHE PATH "")

set(SPHINX_EXECUTABLE /collab/usr/gapps/python/build/spack-toss3.2/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.14-7rci3jkmuht2uiwp433afigveuf4ocnu/bin/sphinx-build CACHE PATH "")
set(DOXYGEN_EXECUTABLE ${GEOSX_TPL_ROOT_DIR}/doxygen/bin/doxygen CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)


set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX /usr/tce/packages/papi/papi-5.4.3 CACHE PATH "")

set(USE_ADDR2LINE ON CACHE BOOL "")

set(ENABLE_OPENMP ON CACHE BOOL "")
set(CUDA_ENABLED OFF CACHE BOOL "")

set(ENABLE_TOTALVIEW_OUTPUT OFF CACHE BOOL "Enables Totalview custom view" FORCE)

set(ENABLE_PETSC OFF CACHE BOOL "Enables PETSc." FORCE)
