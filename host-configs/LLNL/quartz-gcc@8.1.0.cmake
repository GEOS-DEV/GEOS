include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-gcc@8.1.0.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0 CACHE PATH "")

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/libmkl_core.so
                  CACHE STRING "")

# PYGEOSX
set(ENABLE_PYGEOSX ON CACHE BOOL "")
set(Python3_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/linux-rhel7-x86_64-gcc@8.1.0/python@3.8.5/bin/python3 CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)
set(HYPRE_DIR /g/g90/hamon1/geosx/hypre/src/hypre CACHE PATH "" FORCE)
