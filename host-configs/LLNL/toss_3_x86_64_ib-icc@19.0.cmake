
set(CONFIG_NAME "toss_3_x86_64_ib-icc@19.0" CACHE PATH "")

set(COMPILER_BINDIR /usr/tce/packages/intel/intel-19.0.4/compilers_and_libraries_2019.4.227/linux/bin/intel64 )
set(CMAKE_C_COMPILER ${COMPILER_BINDIR}/icc CACHE PATH "")
set(CMAKE_CXX_COMPILER ${COMPILER_BINDIR}/icpc CACHE PATH "")
set(CMAKE_Fortran_COMPILER ${COMPILER_BINDIR}/ifort CACHE PATH "")

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(MPI_HOME             /usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.4 CACHE PATH "")

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/intel/intel-19.0.4/compilers_and_libraries_2019.4.227/linux/mkl)
set(MKL_INCLUDE_DIRS /usr/tce/packages/mkl/mkl-2019.0/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/LLNL/toss_3_x86_64_ib-base.cmake)

set( ENABLE_UNCRUSTIFY OFF CACHE BOOL ""  )