# first Check for RAJA_DIR

if(NOT RAJA_DIR)
    MESSAGE(FATAL_ERROR "Could not find RAJA. RAJA requires explicit RAJA_DIR.")
endif()


#find includes
find_path( RAJA_INCLUDE_DIRS 
           NAMES RAJA/RAJA.hpp
           PATHS  ${RAJA_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( RAJA_LIBRARY 
              NAMES RAJA libRAJA
              PATHS ${RAJA_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set RAJA_FOUND to TRUE
# if all listed variables are TRUE
set(RAJA_FOUND TRUE)

find_package_handle_standard_args(RAJA  DEFAULT_MSG
                                  RAJA_INCLUDE_DIRS
                                  RAJA_LIBRARY )
