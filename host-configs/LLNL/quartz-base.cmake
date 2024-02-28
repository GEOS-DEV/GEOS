###############################################################################
#
# Base configuration for LC Quartz builds
# Calling configuration file must define the following CMAKE variables:
#
# MPI_HOME
#
###############################################################################

# Set GEOS_ROOT_DIR for later use
set(GEOS_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." CACHE PATH "The path to the GEOS root directory" )
message(STATUS "GEOS_ROOT_DIR: ${GEOS_ROOT_DIR}")

# Fortran
set(ENABLE_FORTRAN ON CACHE BOOL "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER ${MPI_HOME}/bin/mpicc CACHE PATH "")
set(MPI_CXX_COMPILER  ${MPI_HOME}/bin/mpicxx CACHE PATH "")
set(MPIEXEC /usr/bin/srun CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")  # test this build 
#set( GEOSX_BUILD_SHARED_LIBS ON CACHE BOOL "" )
#set( GEOSX_BUILD_OBJ_LIBS OFF CACHE BOOL "" )


# PAPI (For TPL caliper builds)
set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX /usr/tce/packages/papi/papi-5.4.3 CACHE PATH "")

# OpenMP
set(ENABLE_OPENMP ON CACHE BOOL "")

# GEOSX specific options
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_PETSC OFF CACHE BOOL "Enables PETSc." FORCE)

# PYGEOSX
set(ENABLE_PYGEOSX ON CACHE BOOL "")
#set(Python3_ROOT_DIR /usr/tce CACHE PATH "")
set(Python3_EXECUTABLE /usr/tce/bin/python3 CACHE PATH "")

# YAPF python formatting
#set(YAPF_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/yapf CACHE PATH "" FORCE)

# Sphinx
#set(SPHINX_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build CACHE PATH "" FORCE)


set(ENABLE_FESAPI OFF CACHE BOOL "" FORCE)

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2022.1.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

# ATS
set(ATS_ARGUMENTS "--machine slurm36"  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
