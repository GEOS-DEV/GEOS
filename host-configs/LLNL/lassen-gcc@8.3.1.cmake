include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-gcc@8.3.1.cmake)

set(GEOSX_TPL_DIR /usr/gapps/GEOSX/thirdPartyLibs/2023-03-15/install-lassen-gcc\@8.3.1-release/ CACHE PATH "" FORCE)
set(ENABLE_CONDUIT ON CACHE BOOL "")

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(FORTRAN_MANGLE_NO_UNDERSCORE OFF CACHE BOOL "")

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1 CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

# Set Inline profiler for MPM to TRUE/FALSE (0/1)
set(MPMPROFILER_OPTION 1) 	# 0 off, 1 something, 2 something else.  doesn't work yet
set(MPM_OPTION 0)      		# 0 default, 1 something, 2 something else.  doesn't work yet

configure_file(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/physicsSolvers/solidMechanics/mpmConfigFiles/config.h.in ${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/physicsSolvers/solidMechanics/temp4P2G/config.h @ONLY)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Build type is Debug")
    message(STATUS "C Flags for Debug configuration: ${CMAKE_C_FLAGS_DEBUG}")
    message(STATUS "CXX Flags for Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

message(STATUS "C Flags for Debug configuration: ${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CXX Flags for Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
