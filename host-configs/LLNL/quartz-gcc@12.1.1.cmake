include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-gcc@12.1.1.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-12.1.1-magic CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)
