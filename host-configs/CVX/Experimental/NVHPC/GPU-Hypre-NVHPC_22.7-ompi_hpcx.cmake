#------------------------------------------------------------------------------
# CUDA support
#------------------------------------------------------------------------------
# _blt_tutorial_surface_cuda_config_start
#
# CVX .cmake for Nvidia nvcc and OpenMPI MPI
#
# file: CPU-GCC_10.2.0-impi_2021.06.cmake
# 
# Michael E. Thomadakis michael.thomadakis@chevron.com

# detect host and name the configuration file
site_name(HOST_NAME)
set(CONFIG_NAME "GPU-Hypre-NVHPC_22.7-ompi_hpcx" CACHE PATH "")
message("## CVX CONFIG_NAME = ${CONFIG_NAME}")
message("## CVX CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}") 

# Default is ON
set(ENABLE_WARNINGS_AS_ERRORS OFF CACHE BOOL "")

# Info
message("## GEOSX_TPL_DIR = ${GEOSX_TPL_DIR}")
message("## GEOSX_DIR = ${GEOSX_DIR}")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "mpif90" CACHE PATH "")
set(MPIEXEC "mpirun" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

set(ENABLE_SPHINX_EXECUTABLE OFF CACHE BOOL "")
set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_ADIAK ON CACHE BOOL "")

set(ENABLE_CUDA ON CACHE BOOL "")
set(ENABLE_CUDA ON CACHE PATH "" FORCE)

## Nvidia HPC TDK 
set(NVHPC_ROOT_DIR "$ENV{NVHPC_ROOT}"  CACHE PATH "")
set(CUDA_TOOLKIT_ROOT_DIR "$ENV{NVHPC_ROOT}/cuda/11.7"  CACHE PATH "")
message("## CUDA_TOOLKIT_ROOT_DIR = ${CUDA_TOOLKIT_ROOT_DIR}")
set(CUDA_TOOLKIT_LIB_DIR "${CUDA_TOOLKIT_ROOT_DIR}/cuda/11.7/lib64:${CUDA_TOOLKIT_ROOT_DIR}/math_libs/lib64"  CACHE PATH "")
message("## CUDA_TOOLKIT_LIB_DIR = ${CUDA_TOOLKIT_LIB_DIR}")
set(LAPACK_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/liblapack_static.a" CACHE STRING "")
message("## LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
set(BLAS_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/libnvblas.so" CACHE STRING "")
message("## BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

# /data/saet/mtml/software/x86_64/nvidia/hpc_sdk/Linux_x86_64/22.7
# CUDA RT libs at 
# /data/saet/mtml/software/aarch64/nvidia/hpc_sdk/Linux_aarch64/22.7/cuda/11.7/lib64/
if(ENABLE_CUDA)
  set(ENABLE_HYPRE_DEVICE CUDA CACHE STRING "")
  message("## ENABLE_CUDA = ${ENABLE_CUDA}; ENABLE_HYPRE_DEVICE = ${ENABLE_HYPRE_DEVICE}")

  # CHECK A100 
  set(CUDA_ARCH sm_80 CACHE STRING "")
  set(CMAKE_CUDA_ARCHITECTURES 80 CACHE STRING "")

  set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
  set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "") 
  set(CMAKE_CUDA_STANDARD 14 CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations " CACHE STRING "") 
#   set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations --library-path=${CUDA_TOOLKIT_LIB_DIR} -L${CUDA_TOOLKIT_LIB_DIR} " CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS_RELEASE " ${CMAKE_CUDA_FLAGS} -O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -fast  " CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0  " CACHE STRING "") 

  set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")

endif()

# CHECK : Most recent version does build. Let's wait for an upgrade on our side.
set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)

if(ENABLE_HYPRE_CUDA)
  set(ENABLE_PETSC OFF CACHE BOOL "") 
  set(ENABLE_TRILINOS OFF CACHE BOOL "") 
  set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "") 
endif()


set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")

# Compilers

# --library-path=${CUDA_TOOLKIT_LIB_DIR}

set(CMAKE_C_COMPILER "nvc" CACHE PATH "")
# set(CMAKE_C_COMPILER "gcc" CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE " -O3 -fast -DNDEBUG -L${CUDA_TOOLKIT_LIB_DIR} " CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")

set(CMAKE_CXX_COMPILER "nvc++" CACHE PATH "")
# set(CMAKE_CXX_COMPILER "g++" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g " CACHE STRING "")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
# set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
# set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

set(CMAKE_Fortran_COMPILER "nvfortran" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE " ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g  " CACHE STRING "")
# set(CMAKE_Fortran_COMPILER "gfortran" CACHE PATH "")
# set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
# set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS " -L${CUDA_TOOLKIT_LIB_DIR} " CACHE STRING "")
# //Flags used by the linker during DEBUG builds.
set(CMAKE_EXE_LINKER_FLAGS_DEBUG " ${CMAKE_EXE_LINKER_FLAGS} " CACHE STRING "")
# //Flags used by the linker during MINSIZEREL builds.
# CMAKE_EXE_LINKER_FLAGS_MINSIZEREL:STRING=

# //Flags used by the linker during RELEASE builds.
set(CMAKE_EXE_LINKER_FLAGS_RELEASE " ${CMAKE_EXE_LINKER_FLAGS}  " CACHE STRING "")

# //Flags used by the linker during RELWITHDEBINFO builds.
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO " ${CMAKE_EXE_LINKER_FLAGS}  " CACHE STRING "")



# enable PAMELA and PVTPackage
set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

# Hypre ON
if(NOT DEFINED ENABLE_HYPRE)
  set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
  set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)
endif()

if(ENABLE_HYPRE)
  set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "" FORCE)
else()
  set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
endif()

if(ENABLE_HYPRE_CUDA)
  set(ENABLE_PETSC OFF CACHE BOOL "") 
  set(ENABLE_TRILINOS OFF CACHE BOOL "") 
  set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "") 
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




# Python options
#set(ENABLE_PYLVARRAY ON CACHE BOOL "")
#set(ENABLE_PYGEOSX ON CACHE BOOL "")
#set(PYTHON_DIR "/share/software/user/open/python/3.6.1" CACHE PATH "")
#set(Python3_EXECUTABLE "/share/software/user/open/python/3.6.1/bin/python3" CACHE PATH "")
