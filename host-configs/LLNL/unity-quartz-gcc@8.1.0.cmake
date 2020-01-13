# CMake for Unity
# cmake executable path: /usr/gapps/GEOSX/thirdPartyLibs/cmake-3.16.2-Linux-x86_64/bin/cmake

include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/LLNL/quartz-gcc@8.1.0.cmake)

# Value "0" combines all sources into one unity file (Default is 8)
set(CMAKE_UNITY_BUILD_BATCH_SIZE 0 CACHE BOOL "")

set(CMAKE_UNITY_BUILD ON CACHE BOOL "")
