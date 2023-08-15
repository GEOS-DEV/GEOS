# CVX .cmake for Intel 1-API compilers and MPI 
#
# file: CPU-GCC_10.2.0-ompi_hpcx.cmake
# 
# Michael E. Thomadakis michael.thomadakis@chevron.com

# detect host and name the configuration file
site_name(HOST_NAME)
set(CONFIG_NAME "CPU-GCC_10.2.0-ompi_hpcx-OMP" CACHE PATH "")

message("## CVX CONFIG_NAME = ${CONFIG_NAME}")
message("## CVX CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}") 

# Default is ON
set(ENABLE_WARNINGS_AS_ERRORS OFF CACHE BOOL "")

 
# # Local mtml settings 
# set(sroot "/data/saet/mtml/src" CACHE PATH "")
# set(iroot "/data/saet/mtml/software/x86_64/RHEL7" CACHE PATH "")
# set(P "GEOSX" CACHE PATH "")
# set(VER "0.2.0" CACHE PATH "")

# # define the path to your compiled installation directory

# if(NOT DEFINED GEOSX_TPL_DIR)
#   set(GEOSX_TPL_ROOT_DIR "${sroot}/${P}/thirdPartyLibs" CACHE PATH "" FORCE)
# else()
#   set(GEOSX_TPL_ROOT_DIR "${iroot}/${P}TPL/${VER}" CACHE PATH "" FORCE)
# endif()
# set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/install-${CONFIG_NAME}-release CACHE PATH "" FORCE)
message("## GEOSX_TPL_DIR = ${GEOSX_TPL_DIR}")

# if(NOT DEFINED GEOSX_DIR)
#   set(GEOSX_ROOT_DIR "${sroot}/${P}/${P}" CACHE PATH "" FORCE)
# else()
#   set(GEOSX_ROOT_DIR "${iroot}/${P}/${VER}" CACHE PATH "" FORCE)
# endif()
# set(GEOSX_DIR ${GEOSX_ROOT_DIR}/install-${CONFIG_NAME}-release CACHE PATH "" FORCE)
message("## GEOSX_DIR = ${GEOSX_DIR}")

# set paths to C, C++, and Fortran compilers. Note that while GEOSX does not contain any Fortran code,
# some of the third-party libraries do contain Fortran code. Thus a Fortran compiler must be specified.
set(CMAKE_C_COMPILER "$ENV{GEOSX_CC}" CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE " -Wno-error -pthread -O3 -Ofast -DNDEBUG " CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

set(CMAKE_CXX_COMPILER "$ENV{GEOSX_CXX}" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

set(CMAKE_Fortran_COMPILER "$ENV{GEOSX_FORT}" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# enable MPI and set paths to compilers and executable.
# Note that the MPI compilers are wrappers around standard serial compilers.
# Therefore, the MPI compilers must wrap the appropriate serial compilers specified
# in CMAKE_C_COMPILER, CMAKE_CXX_COMPILER, and CMAKE_Fortran_COMPILER.
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "$ENV{GEOSX_MPICC}" CACHE PATH "")
set(MPI_CXX_COMPILER "$ENV{GEOSX_MPICXX}" CACHE PATH "")
set(MPI_Fortran_COMPILER "$ENV{GEOSX_MPIFORT}" CACHE PATH "")
set(MPIEXEC "$ENV{GEOSX_MPIRUN}" CACHE PATH "")

# # OpenMP options
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
set(RAJA_ENABLE_OPENMP "ON" CACHE PATH "" FORCE)
# #set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
# #set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")
# disable CUDA and OpenMP
set(CUDA_ENABLED OFF CACHE BOOL "" FORCE)


# # C options
# set(CMAKE_C_COMPILER ${HOST_COMPILER_PATH}/$ENV{GEOSX_CC} CACHE PATH "")
# set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
# set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
# set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

# # C++ options
# set(CMAKE_CXX_COMPILER ${HOST_COMPILER_PATH}/$ENV{GEOSX_CXX} CACHE PATH "")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
# set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
# set(CMAKE_CXX_STANDARD 14 CACHE STRING "")

# # Fortran options
# set(CMAKE_Fortran_COMPILER ${HOST_COMPILER_PATH}/$ENV{GEOSX_FORT} CACHE PATH "")
# set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
# #set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")

# # MPI options
# set(ENABLE_MPI ON CACHE BOOL "")
# set(MPI_ROOT /data_local/sw/spectrum_mpi/10.03.01.00rtm5-rh7_20191114 CACHE PATH "")
# set(MPI_C_COMPILER         ${MPI_ROOT}/bin/$ENV{GEOSX_MPICC}  CACHE PATH "")
# set(MPI_CXX_COMPILER       ${MPI_ROOT}/bin/$ENV{GEOSX_MPICXX} CACHE PATH "")
# set(MPI_Fortran_COMPILER   ${MPI_ROOT}/bin/mpifort CACHE PATH "")
# set(MPIEXEC                ${MPI_ROOT}/bin/$ENV{GEOSX_MPIRUN}  CACHE STRING "")
# set(MPIEXEC_NUMPROC_FLAG   -np CACHE STRING "")
# set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# set(ENABLE_MKL ON CACHE BOOL "")
# set(INTEL_ROOT "$ENV{INTEL_DIR}" )
# set(MKL_ROOT "$ENV{MKLROOT}/mkl" )
# set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
# set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
#                   ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
#                   ${MKL_ROOT}/lib/intel64/libmkl_core.so
#                   ${INTEL_ROOT}/compiler/lib/intel64_lin/libiomp5.so
#                   CACHE STRING "")


# disable CUDA and OpenMP
set(CUDA_ENABLED OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)

# enable PAMELA and PVTPackage
set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")

# Trilinos ON
if(NOT DEFINED ENABLE_TRILINOS)
  set(ENABLE_TRILINOS ON CACHE BOOL "" FORCE)
  # set(ENABLE_TRILINOS "$ENV{ENABLE_TRILINOS}" CACHE BOOL "" FORCE)
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

MESSAGE(STATUS "## GEOSX_LA_INTERFACE = ${GEOSX_LA_INTERFACE}")

# disable Doxygen
set(ENABLE_DOXYGEN OFF CACHE PATH "")

# enable tests
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

# let GEOSX define some third party libraries information for you
message("## CMAKE_CURRENT_LIST_DIR/tpls.cmake = ${CMAKE_CURRENT_LIST_DIR}/tpls.cmake")
include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
message("## CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}")
