include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-gcc@8.3.1.cmake)

# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

# Builds static libraries
set(GEOSX_BUILD_SHARED_LIBS OFF CACHE BOOL "")

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_GEOSX_PTP ON CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")

set(ENABLE_ESSL OFF CACHE BOOL "")

if( ENABLE_ESSL )
    set(ESSL_ROOT_DIR /usr/tcetmp/packages/essl/essl-6.2/ CACHE PATH "" FORCE )
    set(ESSL_INCLUDE_DIRS ${ESSL_ROOT_DIR}/include CACHE PATH "" FORCE)
    set(ESSL_LIB_DIRS ${ESSL_ROOT_DIR}/lib64 CACHE PATH "" FORCE)
    set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.2/lib64/libesslsmpcuda.so
                       /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlsmp.so
                       /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlfmath.so
                       /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlf90_r.so
                       ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so
                       ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
                       ${ESSL_LIB_DIRS}/liblapackforessl.a
                       /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxl.a
                       CACHE PATH "")
else()
    set( BLAS_LIBRARIES /usr/lib64/libblas.so CACHE PATH "")
    set( LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE PATH "")
endif()
set(DOXYGEN_EXECUTABLE /usr/bin/doxygen CACHE PATH "")

set(PETSC_OMP_DIR ${GEOSX_TPL_ROOT_DIR}/omp-links-for-petsc CACHE STRING "")

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "")
