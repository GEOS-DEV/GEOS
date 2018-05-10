##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

#######
# using gcc@4.9.3 compiler spec
#######



set(CONFIG_NAME "quartz-toss_3_x86_64_ib-gcc@7.1.0" CACHE PATH "") 

set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_10_04_13_56_50" CACHE PATH "" )
include("${TPL_DIR}/${CONFIG_NAME}.cmake")

set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-gcc@7.1.0-release" CACHE PATH "")
set(ATK_CMAKE "${ATK_DIR}/lib/cmake" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")

set( GEOSX_TPL_ROOT_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/" CACHE PATH "" )
#set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" FORCE )



#set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
#set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")

#set(ENABLE_SHARED_LIBS ON CACHE BOOL "" FORCE)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON  CACHE BOOL "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )

option( RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" OFF )

option( ENABLE_HYPRE "Enables HYPRE" OFF )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)
