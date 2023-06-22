site_name(HOST_NAME)
set(CONFIG_NAME "cvx-wsl" CACHE PATH "" FORCE)
message( "CONFIG_NAME=${CONFIG_NAME}" )

set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "" FORCE)
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "" FORCE)
set(MPI_Fortran_COMPILER "/usr/bin/mpifort" CACHE PATH "" FORCE)

set(MPIEXEC_EXECUTABLE "/usr/bin/mpirun" CACHE PATH "" FORCE)

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)
set(ENABLE_CALIPER ON CACHE BOOL "")

set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
# If not defined as argument, take from the environment...
if(NOT DEFINED ENABLE_HYPRE)
  set(ENABLE_HYPRE "$ENV{ENABLE_HYPRE}" CACHE BOOL "" FORCE)
endif()
# ... and then check the value.
if(ENABLE_HYPRE)
  set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "" FORCE)
else()
  set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
endif()

# Same pattern
if(NOT DEFINED ENABLE_TRILINOS)
  set(ENABLE_TRILINOS "$ENV{ENABLE_TRILINOS}" CACHE BOOL "" FORCE)
endif()
if(ENABLE_TRILINOS)
  set(GEOSX_LA_INTERFACE "Trilinos" CACHE STRING "" FORCE)
else()
  set(ENABLE_TRILINOS FALSE CACHE BOOL "" FORCE)
endif()

if( (ENABLE_HYPRE AND ENABLE_TRILINOS) OR (NOT ENABLE_TRILINOS AND NOT ENABLE_HYPRE))
  MESSAGE(SEND_ERROR "Exactly one of ENABLE_HYPRE and ENABLE_TRILINOS must be defined.")
  MESSAGE(SEND_ERROR "ENABLE_HYPRE = ${ENABLE_HYPRE}.")
  MESSAGE(SEND_ERROR "ENABLE_TRILINOS = ${ENABLE_TRILINOS}.")
endif()

MESSAGE(STATUS "GEOSX_LA_INTERFACE = ${GEOSX_LA_INTERFACE}")

# set(ENABLE_CUDA "$ENV{ENABLE_CUDA}" CACHE BOOL "" FORCE)
set(ENABLE_CUDA OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
message( "ENABLE_DOXYGEN=${ENABLE_DOXYGEN}" )

if(ENABLE_CUDA)

  set(CMAKE_CUDA_FLAGS "$ENV{CMAKE_CUDA_FLAGS}" CACHE STRING "" FORCE)
  if(NOT CMAKE_CUDA_FLAGS)
    set(CMAKE_CUDA_FLAGS "Unused" CACHE STRING "" FORCE)
  endif()

  set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_TOOLKIT_ROOT_DIR}" CACHE PATH "" FORCE)
  if(NOT CUDA_TOOLKIT_ROOT_DIR)
    set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda" CACHE PATH "" FORCE)
  endif()

  set(CUDA_ARCH "$ENV{CUDA_ARCH}" CACHE STRING "" FORCE)
  if(NOT CUDA_ARCH)
    set(CUDA_ARCH "sm_70" CACHE STRING "" FORCE)
  endif()

  if(ENABLE_HYPRE)
    set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)
  endif()

endif()

# set(GEOSX_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
# set(GEOSX_TPL_DIR "/usr/local/GEOSX/GEOSX_TPL" CACHE PATH "" FORCE )
# if(NOT ( EXISTS "${GEOSX_TPL_DIR}" AND IS_DIRECTORY "${GEOSX_TPL_DIR}" ) )
# Pavel: edit that to provide the path for thirdPartyLibs
#        alternatively (and preferably) use
#          python scripts/config-build.py -hc host-configs/cvx-wsl.cmake -bt Release -D GEOSX_TPL_DIR=/full/path/to/thirdPartyLibs
#set(GEOSX_TPL_DIR "/home/ptls/thirdPartyLibs/install-cvx-wsl-debug" CACHE PATH "" FORCE)
# endif()
include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
