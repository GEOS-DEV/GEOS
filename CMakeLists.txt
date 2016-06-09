cmake_minimum_required(VERSION 3.1)

################################
# ASCTK ToolKit
################################
project(geosx LANGUAGES C CXX)

################################
# Include blt
################################
include(../blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(core/src/cmake/CMakeBasics.cmake)

################################
# Add components
################################
add_subdirectory(core)
add_subdirectory(components)


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
