cmake_minimum_required(VERSION 3.1)

################################
# ASCTK ToolKit
################################
project(geosx LANGUAGES C CXX)

set(CHAI_DIR thirdparty/chai)
################################
# Include blt
################################
include(cmake/blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

################################
# Add components
################################
add_subdirectory(components)


################################
# Add docs
################################
#if (ENABLE_DOCS)
#  add_subdirectory(docs)
#endif()


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
