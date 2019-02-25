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
#module load papi
#export CALI_SERVICES_ENABLE=event:recorder:timestamp:trace:papi
#export CALI_CONFIG_PROFILE=serial-trace
#
#export CALI_TIMER_SNAPSHOT_DURATION=1
#export CALI_TIMER_OFFSET=1
#export CALI_TIMER_TIMESTAMP=1
#export CALI_TIMER_INCLUSIVE_DURATION=1
##################################

##################################
# uses default intel compiler
##################################

set( ENABLE_WARNINGS_AS_ERRORS "OFF" CACHE BOOL "" FORCE)

# c compiler used by spack
set(CMAKE_CXX_COMPILER "CC" CACHE PATH "" FORCE)
set(CMAKE_C_COMPILER "cc" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "ftn" CACHE PATH "" FORCE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -axMIC-AVX512,AVX" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -axMIC-AVX512,AVX" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axMIC-AVX512,AVX" CACHE STRING "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_CXX_COMPILER "${CC}" CACHE PATH "")
set(MPI_C_COMPILER "${cc}" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(GEOSX_TPL_ROOT_DIR "/global/project/projectdirs/m1411/GEOSX/thirdPartyLibs/" CACHE PATH "" )

set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE STRING "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE STRING "" FORCE)

set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "Enables uncrusitfy")

set(HDF5_DIR "/opt/cray/pe/hdf5-parallel/1.10.2.0/cray/8.6" CACHE PATH "" FORCE)

option(RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" On )
set(ENABLE_PAPI "ON" CACHE BOOL "" FORCE)
set(PAPI_PREFIX "/opt/cray/pe/papi/5.5.1.4" CACHE PATH "" FORCE)

option(ENABLE_CHAI "Enables CHAI" ON )
option(ENABLE_RAJA "Enables RAJA" ON )
option(ENABLE_FPARSER "Enables FPARSER" OFF )

set(TRILINOS_TPL_BLAS_LIBRARIES "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/sci_intel_mpi;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mp.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi_mp.a" CACHE PATH "" FORCE )

set(TRILINOS_TPL_BLAS_INCLUDE_DIRS "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )
set(TRILINOS_TPL_BLAS_INCLUDE_DIRS "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )
set(TRILINOS_TPL_LAPACK_LIBRARIES "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mp.a;/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi_mp.a" CACHE PATH "" FORCE)
set(TRILIOS_TPL_LAPACK_INCLUDE_DIRS "/opt/cray/pe/libsci/17.09.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )

set(TPL_CXX_STANDARD "-std=c++14" CACHE STRING "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl  /opt/cray/pe/lib64/libsci_cray.so.5" CACHE STRING "" FORCE)

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

set(ENABLE_PAMELA OFF CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage OFF CACHE BOOL "" FORCE)

#set(CORI_BUILD  CACHE PATH "" FORCE)
set(CUDA_ENABLED      "OFF"       CACHE BOOL "" FORCE)
set(ENABLE_OPENMP     "OFF"       CACHE BOOL "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE STRING "" FORCE)
set(CHAI_ARGS         ""          CACHE STRING "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE BOOL "" FORCE)

# NERSC recommends no optimization flags for intel or cray compilers.
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG" CACHE STRING "" FORCE)

set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)