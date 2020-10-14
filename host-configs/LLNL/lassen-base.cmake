###############################################################################
#
# Base configuration for LC Quartz builds
# Calling configuration file must define the following CMAKE variables:
#
# MPI_HOME
#
###############################################################################

# Fortran
set(ENABLE_FORTRAN OFF CACHE BOOL "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER ${MPI_HOME}/bin/mpicc CACHE PATH "")
set(MPI_CXX_COMPILER ${MPI_HOME}/bin/mpicxx CACHE PATH "")
set(MPIEXEC lrun CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG -n CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# OpenMP
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)

# CUDA
# LvArray sets this to the CMAKE_CXX_COMPILER.
set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")

# ESSL
set(ENABLE_ESSL ON CACHE BOOL "")
set(ESSL_INCLUDE_DIRS /usr/tcetmp/packages/essl/essl-6.2/include CACHE STRING "")
set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.2/lib64/libesslsmpcuda.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlsmp.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlfmath.so
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlf90_r.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
                   ${GEOSX_TPL_ROOT_DIR}/liblapackforesslgeosx.a
                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxl.a
                   CACHE PATH "")

# TPL
set(ENABLE_PAPI OFF CACHE BOOL "")
set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

# GEOSX specific options
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")
set(GEOSX_BUILD_SHARED_LIBS OFF CACHE BOOL "")
set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_GEOSX_PTP ON CACHE BOOL "")
set(ENABLE_PETSC OFF CACHE BOOL "" FORCE )
set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE )
