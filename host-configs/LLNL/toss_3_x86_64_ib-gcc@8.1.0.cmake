set(CONFIG_NAME "toss_3_x86_64_ib-gcc@8.1.0" CACHE PATH "") 

set(CMAKE_C_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/gcc CACHE PATH "")
set(CMAKE_CXX_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/g++ CACHE PATH "")
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-8.1.0/bin/gfortran CACHE PATH "")

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(MPI_HOME             /usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0 CACHE PATH "")

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/libmkl_core.so
                  CACHE STRING "")

set(ENABLE_PETSC ON CACHE BOOL "Enables PETSc." FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/LLNL/toss_3_x86_64_ib-base.cmake)
