include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-clang@14.cmake)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3-clang-14.0.6 CACHE PATH "")

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2022.1.0-gcc-12.1.1/)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

# PYGEOSX
#set(ENABLE_PYGEOSX ON CACHE BOOL "")
#set(Python3_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/linux-rhel7-x86_64-clang@10.0.0/python@3.8.5/bin/python3 CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)
