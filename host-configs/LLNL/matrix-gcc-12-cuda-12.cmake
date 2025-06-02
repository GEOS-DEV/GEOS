include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/matrix-gcc-12-cuda-12.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/matrix-base.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
set(FORTRAN_MANGLE_NO_UNDERSCORE OFF CACHE BOOL "")

# MPI
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

set(ENABLE_CUDA_NVTOOLSEXT OFF CACHE BOOL "")
