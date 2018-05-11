
site_name(HOST_NAME)
set(CONFIG_NAME "${HOST_NAME}-darwin-x86_64-clang@mp4" CACHE PATH "") 
message( "CONFIG_NAME = ${CONFIG_NAME}" )

set(CMAKE_C_COMPILER "/opt/local/bin/clang-mp-5.0" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/opt/local/bin/clang++-mp-5.0" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/opt/local/bin/mpicc-openmpi-clang50" CACHE PATH "")
set(MPI_CXX_COMPILER "/opt/local/bin/mpicxx-openmpi-clang50" CACHE PATH "")
set(MPI_Fortran_COMPILER "/opt/local/bin/mpifort-openmpi-clang50" CACHE PATH "")
set(MPIEXEC "mpirun-openmpi-clang50" CACHE PATH "")

#######################################
# RAJA/CHAI SETUP
#######################################
#set( CHAI_DIR "${CMAKE_SOURCE_DIR}/../../chai" CACHE PATH "")
#set( RAJA_DIR "${CMAKE_SOURCE_DIR}/../../raja" CACHE PATH "")

option( RAJA_ENABLE_TBB "" OFF)
option( ENABLE_CALIPER "Enables CALIPER" OFF )

set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "ON"        CACHE PATH "" FORCE)
