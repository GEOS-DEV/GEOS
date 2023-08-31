include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-icc-19.cmake)

# Fortran
set(CMAKE_Fortran_COMPILER ${COMPILER_DIR}/bin/intel64/ifort CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -march=native -mtune=native -qoverride-limits" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.4 CACHE PATH "")

# GEOSX specific options
set(ENABLE_XML_UPDATES OFF CACHE BOOL "")

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/libmkl_core.so
                  ${COMPILER_DIR}/compiler/lib/intel64/libiomp5.so
                  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)
