###############################################################################
#
# Setup RAVI
# This file defines:
#  RAVI_FOUND - If RAVI was found
#  RAVI_INCLUDE_DIRS - The RAVI include directories

# first Check for RAVI_DIR

if(NOT RAVI_DIR)
    MESSAGE(FATAL_ERROR "Could not find RAVI. RAVI support needs explicit RAVI_DIR")
endif()
set (RAVI_LIBRARY ravi)

include (${RAVI_DIR}/lib/cmake/ravi-targets.cmake)



set (RAVI_INCLUDE_DIRS ${RAVI_DIR}/include)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set RAVI_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(RAVI DEFAULT_MSG
                                  RAVI_INCLUDE_DIRS
                                  RAVI_LIBRARY )

