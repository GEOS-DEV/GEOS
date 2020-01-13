# CMake for Unity
# cmake executable path: /usr/gapps/GEOSX/thirdPartyLibs/cmake-3.16.2-Linux-ppc64le/bin/cmake

include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/LLNL/lassen-clang@upstream.cmake)

# Add missing include path when configuring with CMake 3.16
set(CMAKE_CUDA_FLAGS "-isystem=${CUDA_TOOLKIT_ROOT_DIR}/include ${CMAKE_CUDA_FLAGS}" CACHE STRING "" FORCE)

# Value "0" combines all sources into one unity file (Default is 8)
set(CMAKE_UNITY_BUILD_BATCH_SIZE 0 CACHE BOOL "")

set(CMAKE_UNITY_BUILD ON CACHE BOOL "")
