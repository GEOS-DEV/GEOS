cmake_minimum_required(VERSION 3.1)

################################
# ASCTK ToolKit
################################
project(geosx LANGUAGES C CXX)

################################
# Include blt
################################
include(../asctoolkit/src/blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(src/cmake/CMakeBasics.cmake)

################################
# Add components
################################
add_subdirectory(src)
add_subdirectory(tests)


################################
# Add docs
################################
#if (ENABLE_DOCS)
#  add_subdirectory(docs)
#endif()

################################
# install targets from build
################################s
#install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . OPTIONAL)
#install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . OPTIONAL)
