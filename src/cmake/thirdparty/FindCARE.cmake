###############################################################################
#
# Setup CARE
# This file defines:
#  CARE_FOUND - If CARE was found
#  CARE_INCLUDE_DIRS - The CARE include directories

# first Check for CARE_DIR

if(NOT CARE_DIR)
    MESSAGE(FATAL_ERROR "Could not find CARE. CARE support needs explicit CARE_DIR")
endif()
set (CARE_LIBRARY care)

include (${CARE_DIR}/lib/cmake/care-targets.cmake)



set (CARE_INCLUDE_DIRS ${CARE_DIR}/include)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CARE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CARE DEFAULT_MSG
                                  CARE_INCLUDE_DIRS
                                  CARE_LIBRARY )

