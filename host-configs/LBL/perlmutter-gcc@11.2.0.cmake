include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LBL/perlmutter-gcc@11.2.0.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "" FORCE)

# Fortran
#set(CMAKE_Fortran_COMPILER /opt/cray/pe/craype/2.7.19/bin/ftn CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
#set(FORTRAN_MANGLE_NO_UNDERSCORE OFF CACHE BOOL "")

# MPI
#set(ENABLE_FIND_MPI OFF CACHE BOOL "")
#set(MPI_HOME /opt/cray/pe/mpich/8.1.24/ofi/gnu/9.1 CACHE PATH "")
#set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/perlmutter-base.cmake)
