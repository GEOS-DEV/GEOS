##################################
# cori - KNL build
##################################

##################################
# using default intel compiler
##################################


#set(CONFIG_NAME "quartz-toss_3_x86_64_ib-gcc@7.1.0_noAXOM" CACHE PATH "") 

#set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_10_04_13_56_50" CACHE PATH "" )
#include("${TPL_DIR}/${CONFIG_NAME}.cmake")
#
#set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-gcc@7.1.0-release" CACHE PATH "")
#set(ATK_CMAKE "${ATK_DIR}/lib/cmake" CACHE PATH "")


# c compiler used by spack
#set(SERIAL_CXX_COMPILER "/opt/intel/compilers_and_libraries_2018.0.128/linux/bin/intel64/icpc" CACHE PATH "")
#set(SERIAL_C_COMPILER "/opt/intel/compilers_and_libraries_2018.0.128/linux/bin/intel64/icc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/opt/cray/pe/craype/2.5.12/bin/CC" CACHE PATH "")
set(CMAKE_C_COMPILER "/opt/cray/pe/craype/2.5.12/bin/cc" CACHE PATH "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -std=c++14 -xMIC-AVX512 -DNO_OUTPUT_ON_CORI" CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

#set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-7.1.0" CACHE PATH "")
#set(MPI_HOME             "/opt/intel/compilers_and_libraries_2018.0.128/linux/mpi/intel64" CACHE PATH "")

#set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpiicc"   CACHE PATH "")
#set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpiicpc"  CACHE PATH "")

#set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")


set(MPI_CXX_COMPILER "/opt/cray/pe/craype/2.5.12/bin/CC" CACHE PATH "")
set(MPI_C_COMPILER "/opt/cray/pe/craype/2.5.12/bin/cc" CACHE PATH "")

#set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
#set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
#set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )
#set( GEOSX_TPL_ROOT_DIR "~/Git-Repos/geosx/thirdPartyLibs/" CACHE PATH "" )


set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "Enables uncrusitfy")

#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
#set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )

# We need to use our own chai until chai build issues get resolved
#
set(CHAI_DIR "~/global/homes/v/vargas45/Git-Repos/CHAI/build/install" CACHE PATH "" FORCE )
set(CHAI_INCLUDE_DIRS "~/global/homes/v/vargas45/Git-Repos/CHAI/build/install" CACHE PATH "" FORCE )
set(CHAI_LIBRARY "/global/homes/v/vargas45/Git-Repos/CHAI/build/install/lib/libchai.a" CACHE PATH "" FORCE)

#set(HDF5_DIR "/opt/cray/pe/modulefiles/cray-hdf5-parallel/1.10.0.3" CACHE PATH "" FORCE)
#set(HDF5_DIR "/opt/cray/pe/modulefiles/cray-hdf5/1.8.14" CACHE PATH "" FORCE)
#set(HDF5_DIR "/usr/common/software/modulefiles/hdf5-parallel/1.10.1" CACHE PATH "" FORCE)
 
#set(HDF5_DIR "/global/homes/v/vargas45/Git-Repos/axom/uberenv_libs/spack/opt/spack/linux-x86_64/gcc-7.1.0/hdf5-1.8.16-bptv6mmuhama7dgqbc2zfrye7sk5nirs" CACHE PATH "" FORCE)

#set(TRILINOS_DIR "~/opt/cray/pe/modulefiles/cray-trilinos/12.10.1.1" CACHE PATH "" FORCE)
#set(TRILINOS_DIR "/opt/cray/pe/trilinos/12.10.1.1" CACHE PATH "" FORCE)
#set(TRILINOS_DIR "/opt/cray/pe/trilinos/12.10.1.1/GNU/5.1/x86_64" CACHE PATH "" FORCE)

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


set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl" CACHE PATH "" FORCE)
#set(CORI_BUILD  CACHE PATH "" FORCE)
set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)

#We don't have permission to write to disk... 
#add_compile_options(-D NO_OUTPUT_ON_CORI )
#add_compile_options(NO_OUTPUT_ON_CORI)