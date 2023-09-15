#######################################
#
# Pangea2 - Intel 19 build
#
# SET ENV VARIABLES:
#   export HDF5_USE_FILE_LOCKING=FALSE
#
# Load modules in this order:
#   1) gcc/8.3.0
#   2) intel-compxe/19.1.0.166
#   3) intel-mpi/2020
#   * This is to have the intel compiler use the gcc8 headers instead of the 
#     gcc4.8.5 headers, and to have intel-mpi to wrap intel instead of env cc.
#
########################################


set(CONFIG_NAME "pangea2-icc@19.1.0.166" CACHE PATH "" ) 

set(INTEL_COMPILER_ROOT "/data_local/sw/intel/RHEL7/compilers_and_libraries_2020.0.166/linux" )
set(CMAKE_C_COMPILER "${INTEL_COMPILER_ROOT}/bin/intel64/icc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "${INTEL_COMPILER_ROOT}/bin/intel64/icpc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${INTEL_COMPILER_ROOT}/bin/intel64/ifort" CACHE PATH "")

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(MPI_HOME "${INTEL_COMPILER_ROOT}/mpi/intel64" CACHE PATH "")

#set( BLT_MPI_LINK_FLAGS "-Wl,--enable-new-dtags -Wl,-rpath,/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64/lib/release -Wl,-rpath,/data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mpi/intel64/lib" CACHE PATH "" FORCE )
#set( BLT_MPI_LINK_FLAGS "-Wl,--enable-new-dtags -Wl,-rpath,${MPI_HOME}/lib/release -Wl,-rpath,${MPI_HOME}/lib" CACHE PATH "" FORCE )

set(MKL_ROOT ${INTEL_COMPILER_ROOT}/mkl )
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/TOTAL/pangea2-base.cmake)

