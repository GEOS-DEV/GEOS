include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-clang@9.0.0.cmake)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-4.9.3/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0 CACHE PATH "")

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/libmkl_core.so
                  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)
#set(HYPRE_DIR /usr/workspace/GEOS/GEOSX/geosx_hypre/hypre-mgr-cpu/src/hypre_release CACHE PATH "" FORCE)
