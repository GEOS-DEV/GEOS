# hostconfig for pangea3
# module load cmake/3.14 gcc/8.3 cuda/10.1.243 smpi/10.3.1.0 openblas/0.3.3 lapack/3.9.0

set(CONFIG_NAME "pangea3" CACHE PATH "") 

# Set up the tpls
set(GEOSX_TPL_ROOT_DIR /workrd/SCR/GEOSX/manualClones/thirdPartyLibs CACHE PATH "")
set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/install-${CONFIG_NAME}-release CACHE PATH "")

#set(HOST_COMPILER_PATH /data_local/sw/spack/opt/spack/linux-rhel7-ppc64le/gcc-4.8.5/gcc-8.2.0-suea3az6vbveib7t7d7roevy5y2qeebb/bin )
set(HOST_COMPILER_PATH /data_local/sw/spack/0.13.3/opt/spack/linux-rhel7-power8le/gcc-4.8/gcc-8.3.0-xvcv2vifxgzn2uz52xu3f4njsq5xrc5p/bin )

# C options
set(CMAKE_C_COMPILER ${HOST_COMPILER_PATH}/gcc CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

# C++ options
set(CMAKE_CXX_COMPILER ${HOST_COMPILER_PATH}/g++ CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 14 CACHE STRING "")

# Fortran options
set(CMAKE_Fortran_COMPILER ${HOST_COMPILER_PATH}/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
#set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")








# OpenMP options
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
#set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
#set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

# MPI options
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_ROOT /data_local/sw/spectrum_mpi/10.03.01.00rtm5-rh7_20191114 CACHE PATH "")
set(MPI_C_COMPILER         ${MPI_ROOT}/bin/mpicc  CACHE PATH "")
set(MPI_CXX_COMPILER       ${MPI_ROOT}/bin/mpicxx CACHE PATH "")
set(MPI_Fortran_COMPILER   ${MPI_ROOT}/bin/mpifort CACHE PATH "")
set(MPIEXEC                ${MPI_ROOT}/bin/mpirun  CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG   -np CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# Cuda options
set(ENABLE_CUDA ON CACHE BOOL "")
set(CUDA_TOOLKIT_ROOT_DIR /data_local/sw/cuda/10.1.243 CACHE PATH "")
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "")
set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
set(CUDA_ARCH sm_70 CACHE STRING "")
set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
### The inclusion of -std=c++14 is a workaround for a cuda10/gcc8 bug ###
set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -Xcompiler -std=c++14" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")

# Uncomment this line to make nvcc output register usage for each kernel.
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --resource-usage" CACHE STRING "" FORCE)

# GTEST options
set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)



# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")

set(ENABLE_ESSL OFF CACHE BOOL "")
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
  set( BLAS_LIBRARIES /data_local/sw/openblas/0.3.4/lib/libopenblas.a CACHE PATH "")
  set( LAPACK_LIBRARIES /data_local/sw/lapack/3.9.0/toolchain/rhel-7.6/gcc-8.3/lib64/liblapack.a CACHE PATH "")
endif()

set(ENABLE_DOXYGEN OFF CACHE PATH "")

set(PETSC_OMP_DIR ${GEOSX_TPL_ROOT_DIR}/omp-links-for-petsc CACHE STRING "")

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "")
