###############################################################################
#
# Base configuration for LC Matrix builds
# Calling configuration file must define the following CMAKE variables:
#
# MPI_HOME
#
###############################################################################

set(GEOS_BUILD_OBJ_LIBS OFF CACHE BOOL "")

# Fortran
set(ENABLE_FORTRAN OFF CACHE BOOL "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_HOME "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1" CACHE PATH "")
set(MPI_C_COMPILER ${MPI_HOME}/bin/mpicc CACHE PATH "")
set(MPI_CXX_COMPILER ${MPI_HOME}/bin/mpicxx CACHE PATH "")
set(MPIEXEC srun CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG -n CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# OpenMP
set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)

# CUDA
# LvArray sets this to the CMAKE_CXX_COMPILER.
set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")
set(ENABLE_CUDA_NVTOOLSEXT OFF CACHE BOOL "")

# ESSL
set(ENABLE_ESSL OFF CACHE BOOL "" FORCE )

# TPL
set(ENABLE_PAPI OFF CACHE BOOL "")

# GEOS specific options
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")
set(ENABLE_PETSC OFF CACHE BOOL "" FORCE )
set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE )
set(ENABLE_SUITESPARSE OFF CACHE BOOL "" FORCE )
set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "" FORCE )
set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "" FORCE)
set(ENABLE_HYPRE_MIXINT ON CACHE BOOL "" FORCE)

# Documentation
set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE)
set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)

# YAPF python formatting
#set(YAPF_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python/bin/yapf CACHE PATH "" FORCE)

# PYGEOSX
set(ENABLE_PYGEOSX ON CACHE BOOL "")
set(PYTHON_EXECUTABLE /usr/tce/packages/python/python-3.9.12/bin/python CACHE PATH "")
set(Python3_ROOT_DIR /usr/tce/packages/python/python-3.9.12 CACHE PATH "")
set(Python3_EXECUTABLE /usr/tce/packages/python/python-3.9.12/bin/python3 CACHE PATH "")

# ATS
set(ATS_ARGUMENTS "--ats jsrun_omp --ats jsrun_bind=packed"  CACHE STRING "")
set(USER $ENV{USER} CACHE STRING "")
set(ATS_WORKING_DIR "/p/lustre2/${USER}/integratedTestsGEOS/${CONFIG_NAME}"  CACHE PATH "")
set(ATS_BASELINE_DIR "/p/lustre2/${USER}/integratedTestsGEOS/baselines"  CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
