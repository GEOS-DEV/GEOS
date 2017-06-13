###############################################################################
#
# Setup CHAI
# This file defines:
#  CHAI_FOUND - If CHAI was found
#  CHAI_INCLUDE_DIRS - The CHAI include directories
#  CHAI_LIBRARY - The CHAI library

# first Check for CHAI_DIR

if(NOT CHAI_DIR)
    MESSAGE(FATAL_ERROR "Could not find CHAI. CHAI support needs explicit CHAI_DIR")
endif()

#find includes
find_path( CHAI_INCLUDE_DIRS ManagedArray.hpp
           PATHS  ${CHAI_DIR}/include/chai
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( CHAI_LIBRARY NAMES chai libchai
              PATHS ${CHAI_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CHAI_FOUND to TRUE
# if all listed variables are TRUE
set(CHAI_FOUND TRUE)

find_package_handle_standard_args(CHAI  DEFAULT_MSG
                                  CHAI_INCLUDE_DIRS
                                  CHAI_LIBRARY )
