#######################################
#
# Pangea2 - Intel 19 build
#
# SET ENV VARIABLES:
#   export HDF5_USE_FILE_LOCKING=FALSE
#
# Load modules in this order:
#   1) gcc/8.3.0
#   3) openmpi/2.1.5
#
########################################


set( CONFIG_NAME "pangea2-gcc@8.3.0" CACHE PATH "" ) 

set( CMAKE_C_COMPILER "/data_local/sw/gcc/RHEL7/8.3.0/bin/gcc" CACHE PATH "" )
set( CMAKE_CXX_COMPILER "/data_local/sw/gcc/RHEL7/8.3.0/bin/g++" CACHE PATH "" )
set( CMAKE_Fortran_COMPILER "/data_local/sw/gcc/RHEL7/8.3.0/bin/gfortran" CACHE PATH "" )

set( ENABLE_FORTRAN OFF CACHE BOOL "" FORCE )
set( ENABLE_MPI ON CACHE BOOL "" FORCE )

set( MPI_HOME "/data_local/sw/OpenMPI/RHEL7/2.1.5/gcc/8.3.0" CACHE PATH "" )
set( MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "" )
set( MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "" )
set( MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "" )
set( MPIEXEC              "${MPI_HOME}/bin/mpirun" CACHE PATH "" )
set( MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "" )

set( ENABLE_MKL ON CACHE BOOL "" )
set( MKL_ROOT /data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mkl )
set( MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "" )
set( MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_sequential.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "" )
                  
include( ${CMAKE_CURRENT_LIST_DIR}/../../host-configs/TOTAL/pangea2-base.cmake )
