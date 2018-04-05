
site_name(HOST_NAME)
set(CONFIG_NAME "${HOST_NAME}-darwin-x86_64-clang@mp4" CACHE PATH "") 
message( "CONFIG_NAME = ${CONFIG_NAME}" )


#set(TPL_DIR "${CMAKE_SOURCE_DIR}/../../axom_tpl" CACHE PATH "" )
#message("TPL_DIR = ${TPL_DIR}")
#include("${TPL_DIR}/${CONFIG_NAME}.cmake")
#
#
#
#set(ATK_DIR "${CMAKE_SOURCE_DIR}/../../axom/install-${CONFIG_NAME}-debug" CACHE PATH "")
#set(ATK_CMAKE "${ATK_DIR}/lib/cmake" CACHE PATH "")

set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )


set(CMAKE_C_COMPILER "/opt/local/bin/clang-mp-4.0" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/opt/local/bin/clang++-mp-4.0" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/opt/local/bin/mpicc-openmpi-clang40" CACHE PATH "")
set(MPI_CXX_COMPILER "/opt/local/bin/mpicxx-openmpi-clang40" CACHE PATH "")
set(MPI_Fortran_COMPILER "/opt/local/bin/mpifort-openmpi-clang40" CACHE PATH "")
set(MPIEXEC "mpirun-openmpi-clang40" CACHE PATH "")


include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")


set(GEOSX_LINK_PREPEND_FLAG "-Wl,-force_load" CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "" CACHE PATH "" FORCE)

#set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
#set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
#set( CHAI_DIR "${CMAKE_SOURCE_DIR}/../../chai" CACHE PATH "")
#set( RAJA_DIR "${CMAKE_SOURCE_DIR}/../../raja" CACHE PATH "")
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "OFF"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "OFF"        CACHE PATH "" FORCE)
