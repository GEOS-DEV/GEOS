###############################################################################
#
# Setup TRIBOL WORLD
# This file defines:
#  TRIBOL_WORLD_FOUND - If TRIBOL_WORLD was found
#  TRIBOL_WORLD_INCLUDE_DIRS - The TRIBOL WORLD include directories
#  TRIBOL_WORLD_LIBRARIES - The TRIBOL WORLD libraries

# first Check for TRIBOL_DIR

if(NOT TRIBOL_WORLD_DIR)
    MESSAGE(FATAL_ERROR "Could not find TRIBOL WORLD. TRIBOL WORLD support needs explicit TRIBOL_WORLD_DIR")
endif()

#find includes
find_path( TRIBOL_WORLD_INCLUDE_DIRS tribol/tribol.hpp
           PATHS  ${TRIBOL_WORLD_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)


blt_find_libraries(FOUND_LIBS TRIBOL_WORLD_LIBRARIES 
                   NAMES      slide_world worlds_core rcb_decomp bisectree vista chai tribol
                   PATHS      ${TRIBOL_DIR}/lib )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TRIBOL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TRIBOL_WORLD  DEFAULT_MSG
                                  TRIBOL_WORLD_INCLUDE_DIRS
                                  TRIBOL_WORLD_LIBRARIES )
