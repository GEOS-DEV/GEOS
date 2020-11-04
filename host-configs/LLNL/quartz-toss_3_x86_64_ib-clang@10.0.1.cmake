#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_3_x86_64_ib
# Compiler Spec: clang@10.0.1
# CMake executable path: /usr/tce/packages/cmake/cmake-3.14.5/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-10.0.1/bin/clang" CACHE PATH "")

set(CMAKE_C_FLAGS "-march=native -mtune=native" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-10.0.1/bin/clang++" CACHE PATH "")

set(CMAKE_CXX_FLAGS "-march=native -mtune=native" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.1/bin/mpicxx" CACHE PATH "")

#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(RAJA_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/raja@0.12.1" CACHE PATH "")

set(UMPIRE_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/umpire@4.1.2" CACHE PATH "")

set(CHAI_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/chai@2.2.2" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(HDF5_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/hdf5@1.10.7" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/conduit@master" CACHE PATH "")

set(SILO_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/silo@4.10.2" CACHE PATH "")

set(ADIAK_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/adiak@0.2.1" CACHE PATH "")

set(CALIPER_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/caliper@2.4.0" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/pugixml@1.10" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_MKL ON CACHE BOOL "")

set(MKL_INCLUDE_DIRS "/usr/tce/packages/mkl/mkl-2020.0/include" CACHE PATH "")

set(MKL_LIBRARIES /usr/tce/packages/mkl/mkl-2020.0/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so
                  /usr/tce/packages/mkl/mkl-2020.0/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so
                  /usr/tce/packages/mkl/mkl-2020.0/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64/libmkl_intel_lp64.so
                  /usr/tce/packages/mkl/mkl-2020.0/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64/libmkl_gnu_thread.so
                  /usr/tce/packages/mkl/mkl-2020.0/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64/libmkl_core.so
                  /usr/tce/packages/clang/clang-10.0.1/release/bin/../lib/libomp.so
                  /lib64/libpthread.so
                  /lib64/libm.so
                  /lib64/libdl.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/metis@5.1.0" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/parmetis@4.0.3" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/superlu-dist@6.3.0" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/suite-sparse@5.8.1" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/trilinos@13.0.0" CACHE PATH "")

set(HYPRE_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/hypre@2.20.0" CACHE PATH "")

set(PETSC_DIR "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/petsc@3.14.0" CACHE PATH "")

set(GEOSX_LA_INTERFACE "Trilinos" CACHE PATH "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(ENABLE_PYTHON OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(ENABLE_DOCS OFF CACHE BOOL "")

set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_SPHYNX OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(UNCRUSTIFY_EXECUTABLE "/usr/WS2/corbett5/geosx/uberenv_libs/linux-rhel7-x86_64-clang@10.0.1/uncrustify@0.71.0/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

set(ENABLE_XML_UPDATES OFF CACHE BOOL "")

