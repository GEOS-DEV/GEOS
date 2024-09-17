# file: your-platform.cmake

# detect host and name the configuration file
site_name(HOST_NAME)
set(CONFIG_NAME "ubuntu22" CACHE PATH "")
message("CONFIG_NAME = ${CONFIG_NAME}")

# set paths to C, C++, and Fortran compilers. Note that while GEOS does not contain any Fortran code,
# some of the third-party libraries do contain Fortran code. Thus a Fortran compiler must be specified.
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# enable MPI and set paths to compilers and executable.
# Note that the MPI compilers are wrappers around standard serial compilers.
# Therefore, the MPI compilers must wrap the appropriate serial compilers specified
# in CMAKE_C_COMPILER, CMAKE_CXX_COMPILER, and CMAKE_Fortran_COMPILER.
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/bin/mpifort" CACHE PATH "")
set(MPIEXEC "/usr/bin/mpirun" CACHE PATH "")

# define the path to blas and lapack
#set( BLAS_LIBRARIES /home/rpiazza/lib/lapack-3.11.0-linux-x86_64/libblas.so CACHE PATH "" FORCE )
#set( LAPACK_LIBRARIES /home/rpiazza/lib/lapack-3.11.0-linux-x86_64/liblapack.so  CACHE PATH "" FORCE )

# disable CUDA and OpenMP
set(ENABLE_CUDA OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)

# enable PVTPackage
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

# enable tests
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

# define the path to your compiled installation directory
set(GEOS_TPL_DIR "/home/rpiazza/two-phase/thirdPartyLibs/install-ubuntu22-debug" CACHE PATH "")
#set(GEOSX_TPL_DIR "${GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
# let GEOS define some third party libraries information for you
#include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
