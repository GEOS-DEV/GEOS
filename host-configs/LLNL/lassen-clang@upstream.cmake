include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-clang@upstream.cmake)

#set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/TPLs_boomerAMG-for-elasticity/install-${CONFIG_NAME}-release CACHE PATH "" FORCE )
set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/TPLs_boomerAMG-for-elasticity-cusparse-cublas/install-${CONFIG_NAME}-release CACHE PATH "" FORCE )


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
set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_ESSL ON CACHE BOOL "")
set(ESSL_INCLUDE_DIRS /usr/tcetmp/packages/essl/essl-6.2/include CACHE STRING "")
set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.2/lib64/libesslsmpcuda.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlsmp.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlfmath.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlf90_r.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
                   ${GEOSX_TPL_ROOT_DIR}/liblapackforesslgeosx.a
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxl.a
                   CACHE PATH "")

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "" FORCE )
#set(PETSC_OMP_DIR ${GEOSX_TPL_ROOT_DIR}/omp-links-for-petsc CACHE STRING "")
set(ENABLE_HYPRE ON CACHE BOOL "" FORCE )
