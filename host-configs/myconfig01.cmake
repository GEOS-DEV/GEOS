#----------------------------------------
# GEOSX Quick-Start Configuration
#----------------------------------------
set(CONFIG_NAME "myconfig01") 

#----------------------------------------
# Compilers
#----------------------------------------
# Replace these paths with the correct compiler paths on your system
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "C++ compiler" FORCE)
set(ENABLE_FORTRAN OFF CACHE BOOL "Disable Fortran" FORCE)

#----------------------------------------
# MPI
#----------------------------------------
set(ENABLE_MPI ON CACHE BOOL "Enable MPI" FORCE)
# Replace these paths with your MPI installation
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "MPI C compiler" FORCE)
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "MPI C++ compiler" FORCE)
set(MPIEXEC "/usr/bin/mpirun" CACHE PATH "MPI run command" FORCE)

#----------------------------------------
# BLAS and LAPACK
#----------------------------------------
# Replace these paths with your BLAS/LAPACK library paths
set(BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so" CACHE PATH "BLAS library" FORCE)
set(LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so" CACHE PATH "LAPACK library" FORCE)

#----------------------------------------
# CUDA and OpenMP
#----------------------------------------
set(ENABLE_CUDA OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)

#----------------------------------------
# Third Party Libraries (TPLs)
#----------------------------------------
set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE)
set(ENABLE_CALIPER OFF CACHE BOOL "" FORCE)
set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE)
###
#set(ENABLE_SCOTCH OFF CACHE BOOL "" FORCE)

# Set default TPL directory if not defined
if(NOT (EXISTS "${GEOS_TPL_DIR}" AND IS_DIRECTORY "${GEOS_TPL_DIR}"))
    set(GEOS_TPL_DIR "${CMAKE_SOURCE_DIR}/../../thirdPartyLibs/install-${CONFIG_NAME}-release"
        CACHE PATH "Third-party libraries directory" FORCE)
endif()

# Include additional TPL configuration
include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)

