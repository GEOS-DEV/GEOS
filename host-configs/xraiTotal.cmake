site_name(HOST_NAME)
#set(CONFIG_NAME "${HOST_NAME}-darwin-x86_64-clang@apple-mp" CACHE PATH "") 
#message( "CONFIG_NAME = ${CONFIG_NAME}" )

set(CMAKE_C_COMPILER "/data/PLI/REDHAT7/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/data/PLI/REDHAT7/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/data/PLI/REDHAT7/bin/gfortran" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/data/PLI/REDHAT7/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/data/PLI/REDHAT7/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/data/PLI/REDHAT7/bin/mpifort" CACHE PATH "")
set(MPIEXEC "/data/PLI/REDHAT7/bin/mpirun" CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "OFF"        CACHE PATH "" FORCE)

set(ENABLE_CALIPER "OFF" CACHE PATH "" FORCE )

#option( ENABLE_CALIPER "Enables CALIPER" ON )

set( BLAS_LIBRARIES /data/PLI/REDHAT7/usr/lib64/libblas.so CACHE PATH "" FORCE )
set( LAPACK_LIBRARIES /data/PLI/REDHAT7/usr/lib64/liblapack.so CACHE PATH "" FORCE )

set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)

set(GEOSX_TPL_DIR "/data/PLI/sytuan/GEOSX/thirdPartyLibs/install-xraiTotal_tpls-release" CACHE PATH "" FORCE )

include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
