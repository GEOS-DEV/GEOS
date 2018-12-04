##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

#######
# using gcc@4.9.3 compiler spec
#######



set(CONFIG_NAME "quartz-toss_3_x86_64_ib-gcc@7.1.0_noAXOM" CACHE PATH "") 

#set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_10_04_13_56_50" CACHE PATH "" )
#include("${TPL_DIR}/${CONFIG_NAME}.cmake")
#
#set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-gcc@7.1.0-release" CACHE PATH "")
#set(ATK_CMAKE "${ATK_DIR}/lib/cmake" CACHE PATH "")



# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-7.1.0/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-7.1.0/bin/g++" CACHE PATH "")



set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-7.1.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")



include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")


set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )


set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

set(ENABLE_PAMELA OFF CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage OFF CACHE BOOL "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )

option( RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" OFF )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)
