

set(CONFIG_NAME "quartz-toss_3_x86_64_ib-clang@4.0.0" CACHE PATH "") 

#set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_10_04_13_56_50" CACHE PATH "" )
#include("${TPL_DIR}/${CONFIG_NAME}.cmake")

#set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-clang\@4.0.0-release" CACHE PATH "")

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-4.0.0/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-4.0.0/bin/clang++" CACHE PATH "")


set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)


set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set( GEOSX_TPL_ROOT_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/" CACHE PATH "" )
set(SPHINX_EXECUTABLE "/usr/bin/sphinx-build" CACHE PATH "" FORCE)


#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
#set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )

set(GEOSX_TPL_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/install-toss_3_x86_64_ib-clang@4.0.0-release" CACHE PATH "" FORCE )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "ON"        CACHE PATH "" FORCE)

#######################################
#CALIPER SETUP
#######################################

option( ENABLE_CALIPER "Enables CALIPER" On )
set(GEOSX_TPL_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/install-toss_3_x86_64_ib-clang@4.0.0-release" CACHE PATH "" FORCE )
set(CALIPER_DIR "/g/g17/vargas45/Git-Repo/Caliper/build/" CACHE PATH "" FORCE)