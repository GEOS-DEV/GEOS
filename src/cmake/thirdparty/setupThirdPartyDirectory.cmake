

if(NOT DEFINED GEOSX_TPL_ROOT_DIR)
    get_filename_component(GEOSX_TPL_ROOT_DIR ${CMAKE_SOURCE_DIR}/../../thirdPartyLibs
                           ABSOLUTE CACHE)
endif()

message(STATUS "GEOSX_TPL_ROOT_DIR=${GEOSX_TPL_ROOT_DIR}")

if(NOT DEFINED GEOSX_TPL_DIR)
    get_filename_component(TEMP_DIR "${CMAKE_INSTALL_PREFIX}" NAME)
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        string(REGEX REPLACE "debug" "release" TEMP_DIR2 ${TEMP_DIR})
    endif()
    if(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
        string(REGEX REPLACE "relwithdebinfo" "release" TEMP_DIR2 ${TEMP_DIR})
    endif()
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        set(TEMP_DIR2 ${TEMP_DIR})
    endif()
    set(GEOSX_TPL_DIR "${GEOSX_TPL_ROOT_DIR}/${TEMP_DIR2}")
endif()
message(STATUS "GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")
