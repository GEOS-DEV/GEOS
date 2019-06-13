include("${CMAKE_CURRENT_LIST_DIR}/../src/coreComponents/cxx-utilities/host-configs/ray_blueos_3_ppc64le_ib-clang@upstream.cmake")

# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE)

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE "powerpc64-unknown-linux-gnu" CACHE STRING "" FORCE)

# set(ENABLE_LAPACK_SUITE OFF CACHE BOOL "" FORCE)

# Trilinos options
# set(LAPACK_DIR "/usr/tcetmp/packages/lapack/lapack-3.8.0-P8-gcc-4.9.3/")
# set(TRILINOS_BLAS_LIBRARY_DIRS "${LAPACK_DIR}/lib" CACHE PATH "" FORCE)
# set(TRILINOS_LAPACK_LIBRARY_DIRS "${LAPACK_DIR}/lib" CACHE PATH "" FORCE)
# set(TRILINOS_TPL_BLAS_LIBRARIES "${LAPACK_DIR}/lib/libblas.a" CACHE PATH "" FORCE)
# set(TRILINOS_TPL_BLAS_INCLUDE_DIRS "${LAPACK_DIR}/include/" CACHE PATH "" FORCE)
# set(TRILINOS_TPL_LAPACK_LIBRARIES "${LAPACK_DIR}/lib/liblapack.a" CACHE PATH "" FORCE)
# set(TRILIOS_TPL_LAPACK_INCLUDE_DIRS "${LAPACK_DIR}/include" CACHE PATH "" FORCE)

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
