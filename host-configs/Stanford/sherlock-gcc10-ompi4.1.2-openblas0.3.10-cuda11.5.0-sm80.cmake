set(CONFIG_NAME "sherlock-gcc10-ompi4.1.2-openblas0.3.10-cuda11.5.0-sm80" CACHE PATH "")

set(SOFTWARE_ROOT /share/software/user/open CACHE PATH "")

set(GCC_ROOT "${SOFTWARE_ROOT}/gcc/10.1.0" CACHE PATH "")
set(MPI_ROOT "${SOFTWARE_ROOT}/openmpi/4.1.2" CACHE PATH "")
set(CUDA_HOME "${SOFTWARE_ROOT}/cuda/11.5.0" CACHE PATH "")
set(CMAKE_CUDA_ARCHITECTURES "80" CACHE STRING "")
set(CUDA_ARCH "sm_80" CACHE STRING "")
set(CUDA_TOOLKIT_ROOT_DIR "${CUDA_HOME}" CACHE STRING "")

site_name(HOST_NAME)

# Compilers
set(CMAKE_C_COMPILER       "${GCC_ROOT}/bin/gcc"      CACHE PATH "")
set(CMAKE_CXX_COMPILER     "${GCC_ROOT}/bin/g++"      CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${GCC_ROOT}/bin/gfortran" CACHE PATH "")

# OpenMP options
set(ENABLE_OPENMP OFF CACHE BOOL "")
set(MPI_C_COMPILER       "${MPI_ROOT}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_ROOT}/bin/mpic++"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_ROOT}/bin/mpifort" CACHE PATH "")
set(MPIEXEC              "${MPI_ROOT}/bin/mpirun"  CACHE PATH "")
# MPI options
set(ENABLE_MPI ON CACHE PATH "" FORCE)

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# CUDA options
#set(ENABLE_CUDA OFF)
set(ENABLE_CUDA ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)

message("-- CUDA ARCH ..." ${CUDA_ARCH})
set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")
set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations " CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")

set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")

set(ENABLE_PETSC OFF CACHE BOOL "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "")

# Blas/Lapack options
set(OPENBLAS_ROOT "${SOFTWARE_ROOT}/openblas/0.3.10" CACHE STRING "")
set(BLAS_LIBRARIES "${OPENBLAS_ROOT}/lib/libopenblas.so" CACHE STRING "")
set(LAPACK_LIBRARIES "${OPENBLAS_ROOT}/lib/libopenblas.so" CACHE STRING "")

# Python options
#set(ENABLE_PYLVARRAY ON CACHE BOOL "")
#set(ENABLE_PYGEOSX ON CACHE BOOL "")
#set(PYTHON_DIR "${SOFTWARE_ROOT}/python/3.6.1" CACHE PATH "")
#set(Python3_EXECUTABLE "${SOFTWARE_ROOT}/python/3.6.1/bin/python3" CACHE PATH "")

set(GEOSX_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)


