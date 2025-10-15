set( CONFIG_NAME "ubuntu" ) 

# Set compilers path
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")   # This is typically something like /usr/bin/gcc ... or clang
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "") # This is typically something like /usr/bin/g++ ... or clang++
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# Set paths to mpi
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")    # This is typically something like /usr/bin/mpicc
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "") # This is typically something like /usr/bin/mpicxx
set(MPIEXEC "/usr/bin/mpirun" CACHE PATH "")          # This is typically something like /usr/bin/mpirun

# Set paths to blas and lapack
set( BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so" CACHE PATH "" FORCE )     # This is typically something like /usr/lib64/libblas.so 
set( LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so" CACHE PATH "" FORCE ) # This is typically something like /usr/lib64/liblapack.so

# Cuda and openMP
set( ENABLE_CUDA OFF CACHE PATH "" FORCE )
set( ENABLE_OPENMP OFF CACHE PATH "" FORCE )

# TPLs
set( ENABLE_TRILINOS OFF CACHE PATH "" FORCE )
set( ENABLE_CALIPER OFF CACHE PATH "" FORCE )
set( ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
set( ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )

if(NOT ( EXISTS "${GEOS_TPL_DIR}" AND IS_DIRECTORY "${GEOS_TPL_DIR}" ) )
   set(GEOS_TPL_DIR "${CMAKE_SOURCE_DIR}/../../thirdPartyLibs/install-${CONFIG_NAME}-release" CACHE PATH "" FORCE )
endif()

include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
