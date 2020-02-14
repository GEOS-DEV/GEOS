message( "Entering /src/CMakeLists.txt")

cmake_minimum_required(VERSION 3.9)

################################
# GEOSX
################################
project(geosx LANGUAGES C CXX)
set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")
set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "")
set(CMAKE_ENABLE_EXPORTS ON)

include(cmake/thirdparty/setupThirdPartyDirectory.cmake)



################################
# Include blt
################################
option(ENABLE_BENCHMARKS "Enables benchmarks" ON)
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
set( externalComponentsLinkList "" ) 
add_subdirectory(externalComponents)
add_subdirectory(coreComponents)


################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()


#add_code_check_targets(uncrustify.cfg)


################################
# Deploy headers
################################
install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION .)

################################
# Create header of configuration options
################################
include(cmake/GeosxConfig.cmake)


################################
# Add main
################################
set( extraComponentsLinkList "")
if( ENABLE_PYTHON )
  if( PYTHON_DIR )
    set( extraComponentsLinkList ${extraComponentsLinkList} python_interp )
#    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DGEOSX_USE_PYTHON=1)
    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -L${PYTHON_DIR}/lib)
  else()
    MESSAGE(FATAL_ERROR "ENABLE_PYTHON is true, but PYTHON_DIR is not defined.")
  endif()
else()
#  blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DGEOSX_USE_PYTHON=0)
endif()

if( ENABLE_OPENMP )
  set( extraComponentsLinkList ${extraComponentsLinkList} openmp )
endif()

if( ENABLE_CUDA )
  set( extraComponentsLinkList ${extraComponentsLinkList} cuda )
endif()

blt_add_executable(NAME geosx
                   SOURCES main/main.cpp
                   DEPENDS_ON geosx_core axom
                              ${extraComponentsLinkList}
                              ${externalComponentsLinkList}
                  )

# Removing all third party dependencies from geosx_core target
# Retrieving INTERFACE_LINK_LIBRARIES from geosx_core
get_target_property(linking_list geosx_core INTERFACE_LINK_LIBRARIES)
# Get all CMAKE variables
get_cmake_property(_variable_names VARIABLES)
# Loop over third party libraries
foreach(tpl_library ${thirdPartyLibs})
  string(TOUPPER ${tpl_library} tpl_library_upper)
  set(_tpl_lib "BLT_${tpl_library_upper}_LIBRARIES")
  # Filter to get variables of the form specific BLT targer
  foreach(prop ${_variable_names})
    if(prop MATCHES "^${_tpl_lib}")
      # Save libraries associated to this specific third party dependency
      set(tpl_lib_list ${${prop}})
    endif()
  endforeach()
  # Remove all these items from the main interface linking list
  foreach(_lib ${tpl_lib_list})
    if(${_lib} IN_LIST linking_list)
      list(REMOVE_ITEM linking_list ${_lib})
    endif()
  endforeach()
endforeach()
# Set the new list as INTERFACE_LINK_LIBRARIES variable
set_target_properties(geosx_core PROPERTIES INTERFACE_LINK_LIBRARIES "${linking_list}")

target_include_directories( geosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

set_target_properties( geosx PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

# To change the runtime path during installation
set_target_properties( geosx PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )

install(TARGETS geosx RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

if( ENABLE_MPI AND UNIX AND NOT APPLE AND NOT ENABLE_CUDA )
  add_custom_target( geosx_update_rst_tables
                     ALL
                     COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/updateRstTables.bash ${CMAKE_SOURCE_DIR} > updateRstTables.log
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx
                     )

  add_custom_target( geosx_validate_all_xml_files
                     ALL
                     COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/validateAllXMLFiles.bash ${CMAKE_SOURCE_DIR}
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_update_rst_tables
                     )
endif()

if( ${SPHINX_FOUND} )
  blt_add_sphinx_target( geosx_docs )
endif()

message( "Exiting /src/CMakeLists.txt")
