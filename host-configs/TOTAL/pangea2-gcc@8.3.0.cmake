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

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(MPI_HOME "/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64" CACHE PATH "")

set( BLT_MPI_LINK_FLAGS "-Wl,--enable-new-dtags -Wl,-rpath,/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64/lib/release -Wl,-rpath,/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64/lib" CACHE PATH "" FORCE )

#set(GEOSX_TPL_DIR "/workrd/users/l0505758/geosx/thirdPartyLibs/install-pangea2-gcc@8.3.0-release" CACHE PATH "" FORCE )

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mkl )
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")
                  
include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/TOTAL/pangea2-base.cmake)
