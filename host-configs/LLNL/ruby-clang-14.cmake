include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/ruby-clang-14.cmake)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic CACHE PATH "")

# PYGEOSX
set(ENABLE_PYGEOSX ON CACHE BOOL "")
set(Python3_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python CACHE PATH "")
set(Python3_EXECUTABLE ${Python3_ROOT_DIR}/bin/python3 CACHE PATH "")

# YAPF python formatting
set(YAPF_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/yapf CACHE PATH "" FORCE)

# Sphinx
set(SPHINX_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build CACHE PATH "" FORCE)

# ATS
set(ATS_ARGUMENTS "--machine slurm56"  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/llnl-cpu-base.cmake)
