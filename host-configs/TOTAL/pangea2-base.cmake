#######################################
#
# Pangea2 - Base hostconfig
#
#
########################################


set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpiifort" CACHE PATH "")

set(MPIEXEC              "${MPI_HOME}/bin/mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

#set( GEOSX_TPL_ROOT_DIR /workrd/GEOS/thirdPartyLibs CACHE PATH "")
#set( GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/2019-11-25/install-${CONFIG_NAME}-release CACHE PATH "")

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )
set( ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE )
set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)
set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)
set(ENABLE_MATHPRESSO ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
option( RAJA_ENABLE_TBB "" OFF)
option( ENABLE_CALIPER "Enables CALIPER" ON )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "ON"        CACHE PATH "" FORCE)


