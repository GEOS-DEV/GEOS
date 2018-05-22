message( "Entering /src/CMakeLists.txt")

cmake_minimum_required(VERSION 3.4)

################################
# GEOSX
################################
project(geosx LANGUAGES C CXX)
set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")

################################
# Include blt
################################
include(cmake/blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)
set(CMAKE_DEBUG_POSTFIX "" CACHE STRING "" FORCE)

################################
# Add components
################################
add_subdirectory(components)



################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()


#add_code_check_targets(uncrustify.cfg)


################################
# install targets from build
################################s
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . OPTIONAL)
install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . OPTIONAL)

################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory(thirdparty)

################################
# Create header of configuration options
################################s
include(cmake/GeosxConfig.cmake)

message( "Exiting /src/CMakeLists.txt")
