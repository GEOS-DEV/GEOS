#######################################
#
# Pangea2 - Intel 19 build
#
# SET ENV VARIABLES:
#   export HDF5_USE_FILE_LOCKING=FALSE
#
# Load modules in this order:
#   1) gcc/8.3.0
#   3) intel-mpi/2019U3
#
########################################


set(CONFIG_NAME "pangea2-gcc@8.3.0" CACHE PATH "" ) 

set(CMAKE_C_COMPILER "/data_local/sw/gcc/RHEL7/8.3.0/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/data_local/sw/gcc/RHEL7/8.3.0/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/data_local/sw/gcc/RHEL7/8.3.0/bin/gfortran" CACHE PATH "")


set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpif90" CACHE PATH "")

set(MPIEXEC              "${MPI_HOME}/bin/mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set( BLT_MPI_LINK_FLAGS "-Wl,--enable-new-dtags -Wl,-rpath,/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64/lib/release -Wl,-rpath,/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64/lib" CACHE PATH "" FORCE )


set( ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE )
set(GEOSX_TPL_DIR "/workrd/users/l0505758/geosx/thirdPartyLibs/install-pangea2-gcc@8.3.0-release" CACHE PATH "" FORCE )
#set(SPHINX_EXECUTABLE "/usr/bin/sphinx-build" CACHE PATH "" FORCE)
#set(DOXYGEN_EXECUTABLE "/usr/bin/doxygen" CACHE PATH "" FORCE )

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )


set(ENABLE_PAMELA OFF CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage OFF CACHE BOOL "" FORCE)
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)
set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)
set(ENABLE_MATHPRESSO ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
option( RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" ON )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "ON"        CACHE PATH "" FORCE)


set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mkl )
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")