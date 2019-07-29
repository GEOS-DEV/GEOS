set(CONFIG_NAME "tribol-toss_3_x86_64_ib-gcc@8.1.0" CACHE PATH "") 

set(CMAKE_C_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/gcc CACHE PATH "")
set(CMAKE_CXX_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/g++ CACHE PATH "")
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/gfortran CACHE PATH "")

# set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
# set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "")
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             /usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0 CACHE PATH "")
set(MPI_C_COMPILER       ${MPI_HOME}/bin/mpicc   CACHE PATH "")
set(MPI_CXX_COMPILER     ${MPI_HOME}/bin/mpicxx  CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

set(MPIEXEC              /usr/bin/srun CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )
get_filename_component(ABS_TPL_ROOT_DIR ${GEOSX_TPL_ROOT_DIR} ABSOLUTE)
get_filename_component( INSTALL_DIR_NAME "${CMAKE_INSTALL_PREFIX}" NAME)
set(GEOSX_TPL_DIR ${ABS_TPL_ROOT_DIR}/${INSTALL_DIR_NAME} CACHE PATH "")
set(SPHINX_EXECUTABLE "/usr/bin/sphinx-build" CACHE PATH "" FORCE)
set(DOXYGEN_EXECUTABLE "/usr/bin/doxygen" CACHE PATH "" FORCE )

set(ENABLE_TESTS OFF CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_CALIPER OFF CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_OPENMP OFF CACHE BOOL "")
set(CUDA_ENABLED OFF CACHE BOOL "")

set(ENABLE_MKL OFF CACHE BOOL "")

option(ENABLE_TRIBOL "Enables TRIBOL" ON)
option( RAJA_ENABLE_TBB "" OFF)
option(ENABLE_CHAI "Enables CHAI" ON)

option( ENABLE_CALIPER "Enables CALIPER" OFF )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)

set(RAJA_ENABLE_OPENMP "OFF"        CACHE PATH "" FORCE)
 
set(NUM_PROC "36" CACHE PATH "")
