##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

#######
# using gcc@4.9.3 compiler spec
#######



set(CONFIG_NAME "cab-chaos_5_x86_64_ib-gcc@4.9.3" CACHE PATH "") 

set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_08_25_10_05_28" CACHE PATH "" )
include("${TPL_DIR}/${CONFIG_NAME}.cmake")

set(ATK_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/axom/release" CACHE PATH "")
set(ATK_CMAKE "${ATK_DIR}/lib/cmake" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")




set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

#######################################
# RAJA/CHAI SETUP
#######################################
set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )
option( ENABLE_CALIPER "Enables CALIPER" ON )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)
