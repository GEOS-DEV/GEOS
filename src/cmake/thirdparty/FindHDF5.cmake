###############################################################################
#
# Setup HDF5
# This file defines:
#  HDF5_FOUND - If HDF5 was found
#  HDF5_INCLUDE_DIRS - The HDF5 include directories
#  HDF5_LIBRARY - The HDF5 library

# first Check for HDF5_DIR

if(NOT HDF5_DIR)
    MESSAGE(FATAL_ERROR "Could not find HDF5. HDF5 support needs explicit HDF5_DIR")
endif()

#find includes
find_path( HDF5_INCLUDE_DIRS hdf5.h
           PATHS  ${HDF5_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( HDF5_LIBRARY NAMES hdf5 libhdf5
              PATHS ${HDF5_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set HDF5_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(HDF5  DEFAULT_MSG
                                  HDF5_INCLUDE_DIRS
                                  HDF5_LIBRARY )
