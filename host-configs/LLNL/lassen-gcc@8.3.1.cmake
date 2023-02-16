include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-gcc@8.3.1.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(FORTRAN_MANGLE_NO_UNDERSCORE OFF CACHE BOOL "")

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1 CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/lassen-base.cmake)

# ESSL
set(ENABLE_ESSL ON CACHE BOOL "")
set(ESSL_INCLUDE_DIRS /usr/tcetmp/packages/essl/essl-6.2.1/include CACHE STRING "")
set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.2.1/lib64/libesslsmpcuda.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlsmp.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlfmath.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlf90_r.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
                   /usr/tcetmp/packages/essl/essl-6.2.1/lib64/liblapackforessl.so
                   /usr/tcetmp/packages/essl/essl-6.2.1/lib64/liblapackforessl_.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxl.a
                   CACHE PATH "")