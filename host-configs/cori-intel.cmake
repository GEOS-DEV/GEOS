##################################
# cori - KNL build
#
#Must set the following env variables
# 
#export CRAYPE_LINK_TYPE=dynamic
#export HDF5_USE_FILE_LOCKING=FALSE
#export XTPE_LINK_TYPE=dynamic
#
#And requires CMAKE 3.9 or higher
#export PATH=/global/homes/v/vargas45/Git-Repos/myCMAKE/cmake-3.10.3/bin:$PATH
#
##################################

##################################
# uses default intel compiler
##################################


#set(CONFIG_NAME "quartz-toss_3_x86_64_ib-gcc@7.1.0_noAXOM" CACHE PATH "") 

#set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_10_04_13_56_50" CACHE PATH "" )
#include("${TPL_DIR}/${CONFIG_NAME}.cmake")
#
#set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-gcc@7.1.0-release" CACHE PATH "")
#set(ATK_CMAKE "${ATK_DIR}/lib/cmake" CACHE PATH "")


# c compiler used by spack
set(CMAKE_CXX_COMPILER "CC" CACHE PATH "")
set(CMAKE_C_COMPILER "cc" CACHE PATH "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -std=c++14 -xMIC-AVX512" CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_CXX_COMPILER "${CC}" CACHE PATH "")
set(MPI_C_COMPILER "${cc}" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )

set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "Enables uncrusitfy")

set(HDF5_DIR "/opt/cray/pe/hdf5/1.10.2.0/cray/8.6" CACHE PATH "" FORCE)


option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )

option( RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" OFF )

option( ENABLE_CHAI "Enables CHAI" ON )
option( ENABLE_RAJA "Enables RAJA" ON )
option( ENABLE_FPARSER "Enables FPARSER" OFF )

set(TRILINOS_TPL_BLAS_LIBRARIES "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/sci_intel_mpi;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mp.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi_mp.a" CACHE PATH "" FORCE )

set(TRILINOS_TPL_BLAS_INCLUDE_DIRS "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )
set(TRILINOS_TPL_BLAS_INCLUDE_DIRS "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )
set(TRILINOS_TPL_LAPACK_LIBRARIES "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mp.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi_mp.a" CACHE PATH "" FORCE)
set(TRILIOS_TPL_LAPACK_INCLUDE_DIRS "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )

set(TPL_CXX_STANDARD "-std=c++14" CACHE PATH "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl  /opt/cray/pe/lib64/libsci_cray.so.5" CACHE PATH "" FORCE)

#set(CORI_BUILD  CACHE PATH "" FORCE)
set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)