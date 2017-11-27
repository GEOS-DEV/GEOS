###############################################################################
#
# Setup SILO
# This file defines:
#  SILO_FOUND - If SILO was found
#  SILO_INCLUDE_DIRS - The SILO include directories
#  SILO_LIBRARY - The SILO library

# first Check for SILO_DIR

if(NOT SILO_DIR)
    MESSAGE(FATAL_ERROR "Could not find SILO. SILO support needs explicit SILO_DIR")
endif()

#find includes
find_path( SILO_INCLUDE_DIRS silo.h
           PATHS  ${SILO_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( SILO_LIBRARY NAMES siloh5
              PATHS ${SILO_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SILO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SILO  DEFAULT_MSG
                                  SILO_INCLUDE_DIRS
                                  SILO_LIBRARY )
