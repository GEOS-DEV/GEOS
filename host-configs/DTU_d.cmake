site_name(HOST_NAME)

if(NOT DEFINED CONFIG_NAME)
    set(CONFIG_NAME "${HOST_NAME}-DTU" CACHE STRING "Config name")
endif()

message("CONFIG_NAME = ${CONFIG_NAME}")

set(Python3_EXECUTABLE "/usr/bin/python3" CACHE PATH "")
set(CMAKE_C_COMPILER "gcc" CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-error=array-bounds -Wno-error=dangling-reference -Wno-error=maybe-uninitialized -lgomp -lpthread -lm -ldl" CACHE STRING "")
set(CMAKE_CXX_COMPILER "g++" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-error=array-bounds -Wno-error=dangling-reference -Wno-error=maybe-uninitialized -lgomp -lpthread -lm -ldl" CACHE STRING "")
set(CMAKE_Fortran_COMPILER "gfortran" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-error=array-bounds -Wno-error=dangling-reference -Wno-error=maybe-uninitialized -lgomp -lpthread -lm -ldl" CACHE STRING "")
set(CMAKE_VERBOSE_MAKEFILE ON)

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "mpic++" CACHE PATH "")
set(MPI_Fortran_COMPILER "mpifort" CACHE PATH "")
set(MPIEXEC "mpirun" CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")

set(GEOS_TPL_DIR "/home/behzadh/thirdPartyLibs/install" CACHE PATH "")
set(ENABLE_XML OFF CACHE BOOL "")
set(ENABLE_XML_UPDATES OFF CACHE BOOL "")
set(GEOS_ENABLE_TESTS OFF CACHE BOOL "")
set(ENABLE_EXAMPLES OFF CACHE BOOL "")
set(ENABLE_BENCHMARKS OFF CACHE BOOL "")
set(ENABLE_DOCS OFF CACHE BOOL "")
set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
set(ENABLE_VTK ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

# GTEST
set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
set(gtest_disable_pthreads ON CACHE BOOL "")

# disable most binaries and doc generation
set(ENABLE_TESTS OFF CACHE BOOL "" )
set(DISABLE_UNIT_TESTS ON CACHE BOOL "" )
set(ENABLE_EXAMPLES OFF CACHE BOOL "" )
set(ENABLE_BENCHMARKS OFF CACHE BOOL "" )
set(ENABLE_DOCS OFF CACHE BOOL "" )

set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

set(GEOSX_BUILD_OBJ_LIBS ON CACHE BOOL "")

set(CUDA_ENABLED "OFF" CACHE BOOL "")
set(ENABLE_OPENMP "ON" CACHE BOOL "")

option(ENABLE_CALIPER "Enables CALIPER" ON)

include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT "/opt/intel/oneapi/mkl/latest/")
find_library(MKL_LIBRARIES NAMES mkl_rt PATHS ${MKL_ROOT}/lib/intel64 NO_DEFAULT_PATH)
find_path(MKL_INCLUDE_DIRS NAMES mkl.h PATHS ${MKL_ROOT}/include NO_DEFAULT_PATH)

if(NOT MKL_LIBRARIES OR NOT MKL_INCLUDE_DIRS)
    message(FATAL_ERROR "MKL not found")
endif()

set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

set(CMAKE_SKIP_RULE_DEPENDENCY ON)
set(CMAKE_SKIP_INSTALL_RULES ON)