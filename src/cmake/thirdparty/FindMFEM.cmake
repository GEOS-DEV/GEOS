###############################################################################
#
# Setup MFEM
# This file defines:
#  MFEM_FOUND - If MFEM was found
#  MFEM_INCLUDE_DIRS - The MFEM include directories
#  MFEM_LIBRARY - The MFEM library


# first Check for MFEM_DIR
if(NOT MFEM_DIR)
    MESSAGE(FATAL_ERROR "Could not find MFEM. MFEM support needs explicit MFEM_DIR")
endif()

#find includes
find_path( MFEM_INCLUDE_DIRS mfem.hpp
           PATHS  ${MFEM_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( MFEM_LIBRARY NAMES mfem
              PATHS ${MFEM_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MFEM_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MFEM  DEFAULT_MSG
                                  MFEM_INCLUDE_DIRS
                                  MFEM_LIBRARY )
