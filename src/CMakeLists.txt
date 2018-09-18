message( "Entering /src/CMakeLists.txt")

cmake_minimum_required(VERSION 3.9)

################################
# GEOSX
################################
project(geosx LANGUAGES C CXX)
set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")
set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "")

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
# Add Thirdparty Builtin Libs
################################
add_subdirectory(thirdparty)


################################
# Add components
################################
#add_subdirectory(components)

add_subdirectory(coreComponents)

set( externalComponentsLinkList "" ) 
add_subdirectory(externalComponents)

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
# Create header of configuration options
################################s
include(cmake/GeosxConfig.cmake)









set( extraComponentsListLink "")
if( ENABLE_PYTHON )
  if( PYTHON_DIR )
    set( extraComponentsListLink ${extraComponentsListLink} python_interp )
#    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DGEOSX_USE_PYTHON=1)
    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -L${PYTHON_DIR}/lib)
  else()
    MESSAGE(FATAL_ERROR "ENABLE_PYTHON is true, but PYTHON_DIR is not defined.")
  endif()
else()
#  blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DGEOSX_USE_PYTHON=0)
endif()

if( ENABLE_OPENMP)
  set( extraComponentsListLink ${extraComponentsListLink} openmp )
endif()

blt_add_executable(NAME geosx
                   SOURCES main/main.cpp
                   DEPENDS_ON geosx_core
                              ${extraComponentsLinkList}
                              ${externalComponentsLinkList}
                  )




target_include_directories( geosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


message( "Exiting /src/CMakeLists.txt")
