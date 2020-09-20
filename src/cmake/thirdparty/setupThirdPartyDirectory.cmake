

if( NOT DEFINED GEOSX_TPL_ROOT_DIR )
    get_filename_component( GEOSX_TPL_ROOT_DIR ${CMAKE_SOURCE_DIR}/../../thirdPartyLibs
                           ABSOLUTE  CACHE )
endif()

message(STATUS "GEOSX_TPL_ROOT_DIR=${GEOSX_TPL_ROOT_DIR}")

if( NOT DEFINED GEOSX_TPL_DIR )
    get_filename_component( TEMP_DIR "${CMAKE_INSTALL_PREFIX}" NAME)
    if( CMAKE_BUILD_TYPE MATCHES "Debug" )
        string(REGEX REPLACE "debug" "release" TEMP_DIR2 ${TEMP_DIR})
    endif()
    if( CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo" )
        string(REGEX REPLACE "relwithdebinfo" "release" TEMP_DIR2 ${TEMP_DIR})
    endif()
    if( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(TEMP_DIR2 ${TEMP_DIR})
    endif()
    set( GEOSX_TPL_DIR "${GEOSX_TPL_ROOT_DIR}/${TEMP_DIR2}" )
endif()
message(STATUS "GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")

################################
# uncrustify (This is here instead of SetupGeosxThirdParty.cmake so that BLT will see it)
################################
if ( EXISTS ${UNCRUSTIFY_EXECUTABLE} )
    message( STATUS "Using system uncrustify found at ${UNCRUSTIFY_EXECUTABLE}" )
else()
    set( UNCRUSTIFY_EXECUTABLE "${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify" CACHE PATH "" )
    message( STATUS "Using uncrustify from thirdPartyLibs at ${UNCRUSTIFY_EXECUTABLE}" )
endif()

################################
# Doxygen (This is here instead of SetupGeosxThirdParty.cmake so that BLT will see it)
################################
if ( EXISTS ${DOXYGEN_EXECUTABLE} )
    message( STATUS "Using system Doxygen found at ${DOXYGEN_EXECUTABLE}" )
else()
    set( DOXYGEN_EXECUTABLE "${GEOSX_TPL_DIR}/doxygen/bin/doxygen" CACHE PATH "" )
    message( STATUS "Using Doxygen from thirdPartyLibs at ${DOXYGEN_EXECUTABLE}" )
endif()
