###############################################################################
#
# Setup Vista
# This file defines:
#  VISTA_FOUND - If Vista was found
#  VISTA_INCLUDE_DIRS - The Vista include directories
#  VISTA_LIBRARY - The Vista library

# first Check for VISTA_DIR

if(NOT VISTA_DIR)
    MESSAGE(FATAL_ERROR "Could not find Vista. Vista support needs explicit VISTA_DIR")
endif()

set(VISTA_LIBRARY vista)

include (${VISTA_DIR}/lib/cmake/vista-targets.cmake)

set (VISTA_INCLUDE_DIRS ${VISTA_DIR}/include)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VISTA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VISTA  DEFAULT_MSG
                                  VISTA_INCLUDE_DIRS
                                  VISTA_LIBRARY )
