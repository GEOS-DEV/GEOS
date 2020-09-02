set(CONFIG_NAME "quartz-icc@19.0.4" CACHE PATH "")

set(COMPILER_DIR /usr/tce/packages/intel/intel-19.0.4/compilers_and_libraries_2019.4.227/linux )
set(CMAKE_C_COMPILER ${COMPILER_DIR}/bin/intel64/icc CACHE PATH "")
set(CMAKE_CXX_COMPILER ${COMPILER_DIR}/bin/intel64/icpc CACHE PATH "")
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/bin/intel64/ifort CACHE PATH "")

set(CMAKE_C_FLAGS_DEBUG "-g -O0" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -qoverride-limits" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -qoverride-limits" CACHE STRING "" FORCE)

set(MPI_HOME             /usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.4 CACHE PATH "")

set(ENABLE_XML_UPDATES OFF CACHE BOOL "")

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/libmkl_core.so
                  ${COMPILER_DIR}/compiler/lib/intel64/libiomp5.so
                  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/LLNL/quartz-base.cmake)
