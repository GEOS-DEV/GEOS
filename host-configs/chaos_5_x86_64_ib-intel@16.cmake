##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

#######
# using gcc@4.9.3 compiler spec
#######



set(ATK_ROOT "/usr/gapps/GEOS/asctoolkit" CACHE PATH "")
set(CONFIG_NAME "cab-chaos_5_x86_64_ib-intel@16.0.109" CACHE PATH "") 

set(TPL_DIR "${ATK_ROOT}/thirdparty_libs/2016_10_15_02_06_35" CACHE PATH "" )
set(ATK_DIR "${ATK_ROOT}/install-${CONFIG_NAME}-debug" CACHE PATH "")

message("ATK_DIR=${ATK_DIR}")

include("${CMAKE_CURRENT_LIST_DIR}/hc-defaults.cmake")
include("${TPL_DIR}/${CONFIG_NAME}.cmake")


set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

#######################################
# RAJA/CHAI SETUP
#######################################
set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
set(CALIPER_INSTALL   ""          CACHE PATH "" FORCE)
set(RAJA_ENABLE_TESTS "OFF"       CACHE PATH "" FORCE)
