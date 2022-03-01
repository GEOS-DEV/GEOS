cmake_minimum_required( VERSION 3.9 )

# At the moment we are manually passing the cuda arch flag.
if(${CMAKE_VERSION} VERSION_GREATER "3.18.0") 
  cmake_policy(SET CMP0104 OLD)
endif()

################################
# GEOSX
################################
project( geosx LANGUAGES C CXX )
set( BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard" )
set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "" )
set( CMAKE_ENABLE_EXPORTS ON )

################################
# Include blt
################################
option( ENABLE_BENCHMARKS "Enables benchmarks" ON )
include( cmake/blt/SetupBLT.cmake )
set( BLT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/cmake/blt/ )

# Removes -pthread from BLT CUDA libraries lists
if( ENABLE_CUDA )
  set( _cuda_libraries )
  foreach( _lib ${BLT_CUDA_LIBRARIES} )
    if( "${_lib}" STREQUAL "-pthread" )
      list( APPEND _cuda_libraries "-lpthread" )
    else()
      list( APPEND _cuda_libraries "${_lib}" )
    endif()
  endforeach()

  set( BLT_CUDA_LIBRARIES "${_cuda_libraries}" CACHE STRING "" FORCE )
  set( BLT_CUDA_RUNTIME_LIBRARIES "${_cuda_libraries}" CACHE STRING "" FORCE )
endif()

################################
# Include standard build system logic and options
################################
include( cmake/CMakeBasics.cmake )
set( CMAKE_DEBUG_POSTFIX "" CACHE STRING "" FORCE )


################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory( thirdparty )


################################
# Add components
################################
set( externalComponentsLinkList "" ) 
add_subdirectory( coreComponents )
add_subdirectory( externalComponents )
add_subdirectory( main )


################################
# Add docs
################################
if ( ENABLE_DOCS )
  add_subdirectory( docs )
endif()


################################
# Deploy headers
################################
install( DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . )


################################
# Generate version information
################################
include( cmake/GeosxVersion.cmake )
message( STATUS "Configuring GEOSX version ${GEOSX_VERSION_FULL}" )


################################
# Create header of configuration options
################################
include( cmake/GeosxConfig.cmake )


################################
# Add main
################################
set( extraComponentsLinkList "" )

if( ENABLE_OPENMP )
  set( extraComponentsLinkList ${extraComponentsLinkList} openmp )
endif()

if( ENABLE_CUDA )
  set( extraComponentsLinkList ${extraComponentsLinkList} cuda )
endif()

if ( GEOSX_BUILD_SHARED_LIBS )
  set( extraComponentsLinkList ${extraComponentsLinkList} geosx_core )
else()
  set( extraComponentsLinkList ${extraComponentsLinkList} ${geosx_core_libs} )
endif()

blt_add_executable( NAME geosx
                    SOURCES main/main.cpp
                    DEPENDS_ON ${extraComponentsLinkList}
                               ${externalComponentsLinkList}
                  )

# Seems to be required on some CMake versions (e.g. 3.16) to get enforce device linking
if( ENABLE_HYPRE_CUDA )
  set_target_properties( geosx PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS TRUE )
endif()

# Removing all transitive link dependencies from geosx_core target to circumvent
# the BLT behavior which imposes all dependencies to be public
#set_target_properties(geosx_core PROPERTIES INTERFACE_LINK_LIBRARIES "")

target_include_directories( geosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

# To change the runtime path during installation
set_target_properties( geosx PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
set_target_properties( geosx PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

install(TARGETS geosx RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

if( ENABLE_XML_UPDATES AND ENABLE_MPI AND UNIX AND NOT CMAKE_HOST_APPLE AND NOT ENABLE_CUDA )

  set(SCHEMA_DIR ${CMAKE_SOURCE_DIR}/coreComponents/schema)
  set(SCRIPT_DIR ${CMAKE_SOURCE_DIR}/../scripts)

  add_custom_target( geosx_generate_schema
                     ALL
                     COMMAND geosx -s ${SCHEMA_DIR}/schema.xsd >generate_schema.log 2>&1 || (cat generate_schema.log && exit 1)
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx
                     COMMENT "Generating XML schema"
                   )

  add_custom_target( geosx_update_rst_tables
                     ALL
                     COMMAND python ${SCRIPT_DIR}/SchemaToRSTDocumentation.py ${SCHEMA_DIR} >update_rst_tables.log 2>&1 || (cat update_rst_tables.log && exit 1)
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_generate_schema
                     COMMENT "Generating schema-dependent RST files"
                   )

  add_custom_target( geosx_validate_all_xml_files
                     ALL
                     COMMAND bash ${SCRIPT_DIR}/validateXMLFiles.bash -g ${SCHEMA_DIR}/schema.xsd ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../examples ${CMAKE_SOURCE_DIR}/../inputFiles
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_update_rst_tables
                     COMMENT "Validating all XML files in the repository against the schema"
                   )
endif()

if( SPHINX_FOUND )
  blt_add_sphinx_target( geosx_docs )
endif()


install( DIRECTORY ../examples/ 
         DESTINATION examples )

install( DIRECTORY ../scripts/ 
         DESTINATION scripts )

if( ENABLE_PYGEOSX )
  add_subdirectory( pygeosx )
endif()
