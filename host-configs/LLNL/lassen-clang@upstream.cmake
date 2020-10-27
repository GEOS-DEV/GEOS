include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-clang@upstream.cmake)

set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/TPLs_boomerAMG-for-elasticity-2020-10-22/install-${CONFIG_NAME}-release CACHE PATH "" FORCE )
#set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/TPLs_boomerAMG-for-elasticity-hypredebug/install-${CONFIG_NAME}-release CACHE PATH "" FORCE )
#set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/TPLs_boomerAMG-for-elasticity-cusparse-cublas/install-${CONFIG_NAME}-release CACHE PATH "" FORCE )


# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/xl/xl-beta-2019.06.20/bin/xlf_r CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -qarch=pwr9 -qtune=pwr9" CACHE STRING "")
set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")
set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-upstream-2019.03.26 CACHE PATH "")
set(MPI_Fortran_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-beta-2019.06.20/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/lassen-base.cmake)
