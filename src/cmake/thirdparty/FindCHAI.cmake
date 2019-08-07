###############################################################################
#
# Setup CHAI
# This file defines:
#  CHAI_FOUND - If CHAI was found
#  CHAI_INCLUDE_DIRS - The CHAI include directories
#  CHAI_LIBRARY - The CHAI library

# first Check for CHAI_DIR

if (NOT EXISTS ${CHAI_DIR})
    message(STATUS "Using chai from thirdPartyLibs")
    set(CHAI_DIR "${GEOSX_TPL_DIR}/chai" CACHE PATH "")
endif()

if (EXISTS ${CHAI_DIR})
    message(STATUS "Using chai found at ${CHAI_DIR}")
    
    find_library( CHAI_LIBRARY NAMES chai libchai libumpire libumpire_op libumpire_resource libumpire_strategy libumpire_tpl_judy libumpire_util
                  PATHS ${CHAI_DIR}/lib
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH)

    include("${CHAI_DIR}/share/umpire/cmake/umpire-targets.cmake")
    include("${CHAI_DIR}/share/chai/cmake/chai-targets.cmake")

    # handle the QUIETLY and REQUIRED arguments and set CHAI_FOUND to TRUE
    # if all listed variables are TRUE
    set(CHAI_INCLUDE_DIRS ${CHAI_DIR}/include)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args( CHAI  DEFAULT_MSG
                                        CHAI_INCLUDE_DIRS
                                        CHAI_LIBRARY )

    if (NOT CHAI_FOUND)
        message(FATAL_ERROR "CHAI not found in ${CHAI_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME chai
                          INCLUDES ${CHAI_INCLUDE_DIRS}
                          LIBRARIES ${CHAI_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON )

    set(ENABLE_CHAI ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} chai)
else()
    set(CHAI_FOUND FALSE)
    set(ENABLE_CHAI OFF CACHE BOOL "")
    message(STATUS "Not using chai")
endif()
