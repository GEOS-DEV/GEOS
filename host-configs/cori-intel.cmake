##################################
# cori - KNL build
#
# Must set the following env variables
# 
# export CRAYPE_LINK_TYPE=dynamic
# export HDF5_USE_FILE_LOCKING=FALSE
# export XTPE_LINK_TYPE=dynamic
#
# And load the following
# module load cmake/3.11.4
# module load gcc/7.3.0
# module swap intel/18.0.1.163 intel/19.0.3.199
# module swap craype/2.5.15 craype/2.5.18
# module swap cray-libsci/18.07.1 cray-libsci/19.02.1
##################################

set( ENABLE_WARNINGS_AS_ERRORS "OFF" CACHE BOOL "" FORCE)

set(CMAKE_CXX_COMPILER "CC" CACHE PATH "" FORCE)
set(CMAKE_C_COMPILER "cc" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "ftn" CACHE PATH "" FORCE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -axMIC-AVX512,AVX" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -axMIC-AVX512,AVX" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axMIC-AVX512,AVX" CACHE STRING "" FORCE)

# NERSC recommends no optimization flags for intel or cray compilers.
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG" CACHE STRING "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "" FORCE)
set(MPI_CXX_COMPILER "${CC}" CACHE PATH "")
set(MPI_C_COMPILER "${cc}" CACHE PATH "")
set(MPI_Fortran_COMPILER "ftn" CACHE PATH "" FORCE)
set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(GEOSX_TPL_DIR "/global/homes/c/corbett5/GEOSX/tpls/install-cori-intel-release" CACHE PATH "" )

set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE STRING "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE STRING "" FORCE)

set(ENABLE_SPHINX_EXECUTABLE OFF CACHE BOOL "")
set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_PAMELA OFF CACHE BOOL "")
set(ENABLE_PVTPackage OFF CACHE BOOL "")
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX "/opt/cray/pe/papi/5.7.0.1" CACHE PATH "" FORCE)

set(ENABLE_OPENMP ON CACHE BOOL "")
set(CUDA_ENABLED OFF CACHE BOOL "")

set(ENABLE_TOTALVIEW_OUTPUT OFF CACHE BOOL "Enables Totalview custom view" FORCE)

# set(TRILINOS_TPL_BLAS_LIBRARIES "/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel.a;/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/sci_intel_mpi;/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel_mp.a;/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi_mp.a" CACHE PATH "" FORCE )
# set(TRILINOS_TPL_BLAS_INCLUDE_DIRS "/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )
# set(TRILINOS_TPL_LAPACK_LIBRARIES "/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel.a;/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi.a;/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel_mp.a;/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/lib/libsci_intel_mpi_mp.a" CACHE PATH "" FORCE)
# set(TRILIOS_TPL_LAPACK_INCLUDE_DIRS "/opt/cray/pe/libsci/19.02.1/INTEL/16.0/x86_64/include" CACHE PATH "" FORCE )
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl /opt/cray/pe/lib64/libsci_cray.so.5" CACHE STRING "" FORCE)

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

set(ENABLE_PETSC OFF CACHE BOOL "")
