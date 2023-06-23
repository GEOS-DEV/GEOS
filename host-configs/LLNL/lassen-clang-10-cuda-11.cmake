include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-clang-10.cmake)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/xl/xl-2022.08.19-cuda-11.8.0/bin/xlf_r  CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -qarch=pwr9 -qtune=pwr9" CACHE STRING "")
set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")
set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-10.0.1-gcc-8.3.1 CACHE PATH "")
set(MPI_Fortran_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19-cuda-11.8.0/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/lassen-base.cmake)

set(ENABLE_CUDA_NVTOOLSEXT OFF CACHE BOOL "")

set(ENABLE_ESSL ON CACHE BOOL "" FORCE )
set(ESSL_INCLUDE_DIRS /usr/tcetmp/packages/essl/essl-6.3.0.2/include CACHE STRING "" FORCE )
set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/libesslsmpcuda.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlsmp.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlfmath.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlf90_r.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublasLt.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
                   /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/liblapackforessl.so
                   /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/liblapackforessl_.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxl.a
                   CACHE PATH "" FORCE )