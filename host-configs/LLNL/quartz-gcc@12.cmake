include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-gcc@12.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-12.1.1 CACHE PATH "")

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2022.1.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")


include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)