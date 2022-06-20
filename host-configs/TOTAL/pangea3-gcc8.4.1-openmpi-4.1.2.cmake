# hostconfig for pangea3
#
# export MODULEPATH=/workrd/SCR/NUM/geosx_num/module_files:$MODULEPATH
# module load cmake/3.21.4 gcc/8.4.1 cuda/11.0.3 ompi/4.1.2 openblas/0.3.18 python4geosx/p3/gcc8.4.1-ompi4.1.2
#
set(CONFIG_NAME "pangea3-gcc8.4.1-ompi-4.1.2" CACHE PATH "") 

# Set up the tpls
if (NOT DEFINED GEOSX_TPL_DIR)
  message(FATAL_ERROR "You must set GEOSX_TPL_DIR with -D GEOSX_TPL_DIR=")
endif ()

# C options
set(CMAKE_C_COMPILER gcc CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

# C++ options
set(CMAKE_CXX_COMPILER g++ CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 14 CACHE STRING "")

# Fortran options
set(CMAKE_Fortran_COMPILER gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
#set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")


# OpenMP options
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
#set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
#set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

# MPI options
if (DEFINED ENV{MPI_ROOT})
  set(ENABLE_MPI ON CACHE BOOL "")
  set(MPI_C_COMPILER         $ENV{MPI_ROOT}/bin/mpicc  CACHE PATH "")
  set(MPI_CXX_COMPILER       $ENV{MPI_ROOT}/bin/mpicxx CACHE PATH "")
  set(MPI_Fortran_COMPILER   $ENV{MPI_ROOT}/bin/mpifort CACHE PATH "")
  set(MPIEXEC                $ENV{MPI_ROOT}/bin/mpirun  CACHE STRING "")
  set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")
else()
  message(FATAL_ERROR "You must have MPI_ROOT variable set, we advise loading module ompi/4.1.2")
endif()

# Cuda options
if (DEFINED ENV{CUDA_ROOT})
  set(ENABLE_CUDA ON CACHE BOOL "")
  set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_ROOT} CACHE PATH "")
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "")
  set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
  set(CUDA_ARCH sm_70 CACHE STRING "")
  set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "")
  set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
  ### The inclusion of -std=c++14 is a workaround for a cuda10/gcc8 bug ###
  set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -Xcompiler -std=c++14" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")

  # Uncomment this line to make nvcc output register usage for each kernel.
  # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --resource-usage" CACHE STRING "" FORCE)
else()
  message(FATAL_ERROR "You must have CUDA_ROOT environment variable set, we advise loading module cuda/11.0.3")
endif()

# GTEST options
set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)



# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

set(GEOSX_BUILD_SHARED_LIBS OFF CACHE BOOL "")
set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")

set(ENABLE_ESSL OFF CACHE BOOL "")

set(ENABLE_PYGEOSX ON CACHE BOOL "")

if( ENABLE_ESSL )
    set(ESSL_DIR /data_local/sw/essl/6.1 CACHE PATH "")
    set(XL_DIR /data_local/sw/xl/16.1.1.3-190404 CACHE PATH "" )
    set(ESSL_INCLUDE_DIRS ${ESSL_DIR}/include CACHE PATH "")   
    set(ESSL_LIBRARIES ${ESSL_DIR}/lib64/libesslsmpcuda.so
    		       ${XL_DIR}/xlsmp/5.1.1/lib/libxlsmp.so	
		       ${XL_DIR}/xlf/16.1.1/lib/libxlfmath.so
       		       ${XL_DIR}/xlf/16.1.1/lib/libxlf90_r.so
                       ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas_static.a
                       ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart_static.a
                       ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcufft_static.a
                       ${GEOSX_TPL_ROOT_DIR}/liblapackforessl.a
		       ${XL_DIR}/xlC/16.1.1/lib/libxl.a
                       CACHE PATH "")
else()
  set(BLAS_LIBRARIES /data_local/sw/spack/0.17.0/opt/spack/linux-rhel8-power9le/gcc-8.4.1/openblas-0.3.18-udwdz2i4a3zcoyjl63h2wlsoacmifvwk/lib/libopenblas.a)
  set(LAPACK_LIBRARIES /data_local/sw/spack/0.17.0/opt/spack/linux-rhel8-power9le/gcc-8.4.1/openblas-0.3.18-udwdz2i4a3zcoyjl63h2wlsoacmifvwk/lib/libopenblas.a)
endif()

set(ENABLE_DOXYGEN OFF CACHE PATH "")

set(PETSC_OMP_DIR ${GEOSX_TPL_ROOT_DIR}/omp-links-for-petsc CACHE STRING "")

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "")

include( ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake )
