#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_ppc64le_ib
# Compiler Spec: clang@=14.0.6
# CMake executable path: /usr/tce/packages/cmake/cmake-3.23.1/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-14.0.6/bin/clang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-14.0.6/bin/clang++" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2023.06.28-cuda-11.8.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2023.06.28-cuda-11.8.0/bin/mpicxx" CACHE PATH "")

#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA ON CACHE BOOL "")

set(CMAKE_CUDA_STANDARD "17" CACHE PATH "")

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-11.8.0/" CACHE PATH "")

set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE PATH "")

set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")

set(CMAKE_CUDA_ARCHITECTURES "cuda_arch=70" CACHE STRING "")

set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -arch cuda_arch=70" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -Xcompiler -O3 -DNDEBUG" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo -O3 -Xcompiler -O3 -DNDEBUG" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0 " CACHE STRING "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(RAJA_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/raja-develop" CACHE PATH "")

set(UMPIRE_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/umpire-develop" CACHE PATH "")

set(CHAI_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/chai-develop" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(HDF5_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/hdf5-1.14.2" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/conduit-0.8.8" CACHE PATH "")

set(ENABLE_SILO OFF CACHE BOOL "")

set(ENABLE_ADIAK OFF CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/caliper-2.10.0" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/pugixml-1.13" CACHE PATH "")

set(ENABLE_VTK OFF CACHE BOOL "")

set(FMT_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/fmt-10.1.0" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_ESSL ON CACHE BOOL "")

set(ESSL_INCLUDE_DIRS "/usr/tcetmp/packages/essl/essl-6.3.0.2/include" CACHE PATH "")

set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/libessl.so
                   /usr/tce/packages/cuda/cuda-11.8.0/lib64/libcudart.so CACHE STRING "")

set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/metis-5.1.0" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/parmetis-4.0.3" CACHE PATH "")

set(ENABLE_SCOTCH OFF CACHE BOOL "")

set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "")

set(ENABLE_SUITESPARSE OFF CACHE BOOL "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/tobin6/geosx/develop/scripts/uberenv/uberenv_libs/spack/opt/spack/install-linux-rhel8-ppc64le-clang-14.0.6/hypre-develop" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")

set(ENABLE_PETSC OFF CACHE BOOL "")

set(GEOSX_LA_INTERFACE "Hypre" CACHE PATH "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(ENABLE_PYGEOSX OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(ENABLE_DOCS OFF CACHE BOOL "")

set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_SPHINX OFF CACHE BOOL "")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

set(ENABLE_XML_UPDATES OFF CACHE BOOL "")

