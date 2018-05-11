

set(CONFIG_NAME "quartz-toss_3_x86_64_ib-clang@4.0.0" CACHE PATH "") 

set(TPL_DIR "/usr/gapps/GEOS/geosx/2017_10_04_13_56_50" CACHE PATH "" )
include("${TPL_DIR}/${CONFIG_NAME}.cmake")

set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-clang\@4.0.0-release" CACHE PATH "")


set( GEOSX_TPL_ROOT_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/" CACHE PATH "" FORCE)



#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
#set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )

option( RAJA_ENABLE_TBB "" OFF)

option( ENABLE_CALIPER "Enables CALIPER" OFF )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "ON"        CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "ON"        CACHE PATH "" FORCE)
