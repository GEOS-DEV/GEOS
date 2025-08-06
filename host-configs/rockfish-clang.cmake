set( CONFIG_NAME "rockfish" ) 

# Set compilers path
set(CMAKE_C_COMPILER "/data/apps/extern/spack_on/gcc/9.3.0/llvm/17.0.4-y742actqtds5jan3z6a3yu4fuyp7edtw/bin/clang" CACHE PATH "")   # This is typically something like /usr/bin/gcc ... or clang
set(CMAKE_CXX_COMPILER "/data/apps/extern/spack_on/gcc/9.3.0/llvm/17.0.4-y742actqtds5jan3z6a3yu4fuyp7edtw/bin/clang++" CACHE PATH "") # This is typically something like /usr/bin/g++ ... or clang++
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# Set paths to mpi
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/openmpi-3.1.6-rk3nyoehbq3pke4zy4hn7unns3ujtutx/bin/mpicc" CACHE PATH "")    # This is typically something like /usr/bin/mpicc
set(MPI_CXX_COMPILER "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/openmpi-3.1.6-rk3nyoehbq3pke4zy4hn7unns3ujtutx/bin/mpicxx" CACHE PATH "") # This is typically something like /usr/bin/mpicxx
set(MPIEXEC "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/openmpi-3.1.6-rk3nyoehbq3pke4zy4hn7unns3ujtutx/bin/mpirun" CACHE PATH "")          # This is typically something like /usr/bin/mpirun

# Set paths to blas and lapack
set( BLAS_LIBRARIES "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/openblas-0.3.10-hnsda7rer5flvldg5upobbsnfiw437jp/lib/libopenblas.so" CACHE PATH "" FORCE )     # This is typically something like /usr/lib64/libblas.so 
set( LAPACK_LIBRARIES "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/netlib-lapack-3.8.0-6a22x27vpnialeqoxzyhzhaeaqhvtbku/lib/liblapack.so" CACHE PATH "" FORCE ) # This is typically something like /usr/lib64/liblapack.so

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
