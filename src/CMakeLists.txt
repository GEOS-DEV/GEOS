
cmake_minimum_required( VERSION 3.23 )

# At the moment we are manually passing the cuda arch flag.
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.23.0")
  cmake_policy(SET CMP0104 OLD) # when using nvcc populate CMAKE_CUDA_ARCHITECTURES, raise error if we can't
endif()
cmake_policy(SET CMP0074 NEW) # dont ignore <PackageName>_ROOT env vars when searching for packages via find_package()
cmake_policy(SET CMP0066 NEW) # use CMAKE_<LANG>_FLAGS_<BUILD_TYPE> for try_compile() instead of only CMAKE_<LANG>_FLAGS
cmake_policy(SET CMP0056 NEW) # use CMAKE_EXE_LINKER_FLAGS in try_compile() in addition to CMAKE_<LANG>_FLAGS

################################
# GEOSX
################################
project( geosx LANGUAGES C CXX )
include(GNUInstallDirs)
set( BLT_CXX_STD "c++17" CACHE STRING "Version of C++ standard" )
if(CMAKE_BUILD_TYPE EQUAL "Debug")
  set( ENABLE_WARNINGS_AS_ERRORS "OFF" CACHE PATH "")
else()
  set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "")
endif()
set( CMAKE_ENABLE_EXPORTS ON )

set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CMAKE_COMMAND} -E time")

################################
# Include blt
################################
option( ENABLE_BENCHMARKS "Enables benchmarks" ON )
include( cmake/blt/SetupBLT.cmake )
set( BLT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/cmake/blt/ )

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

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

if( ENABLE_CUDA_NVTOOLSEXT )
  set( extraComponentsLinkList ${extraComponentsLinkList} CUDA::nvToolsExt )
endif()

if( ENABLE_HIP )
  set( extraComponentsLinkList ${extraComponentsLinkList} blt::hip )
endif()

if ( GEOSX_BUILD_SHARED_LIBS )
  set( extraComponentsLinkList ${extraComponentsLinkList} geosx_core )
else()
  set( extraComponentsLinkList ${extraComponentsLinkList} ${geosx_core_libs} )
endif()

blt_add_executable( NAME geosx
                    SOURCES main/main.cpp
                    DEPENDS_ON ${extraComponentsLinkList}
                               ${externalComponentsLinkList} )

# Seems to be required on some CMake versions (e.g. 3.16) to get enforce device linking
if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA"  )
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
# Duplicate the installed `geosx` executable and rename it `geos` to let the end users adapt their scripts.
install( FILES ${CMAKE_INSTALL_PREFIX}/bin/geosx
         DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
         RENAME geos
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

if( ENABLE_XML_UPDATES AND ENABLE_MPI AND UNIX AND NOT CMAKE_HOST_APPLE AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )

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
                     COMMAND python3 ${SCRIPT_DIR}/SchemaToRSTDocumentation.py ${SCHEMA_DIR} >update_rst_tables.log 2>&1 || (cat update_rst_tables.log && exit 1)
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

if( ENABLE_PYGEOSX )
   add_subdirectory( pygeosx )
endif()

if( SPHINX_FOUND )
  blt_add_sphinx_target( geosx_docs )
endif()


install( DIRECTORY ../examples/
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/examples )

install( DIRECTORY ../scripts/
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/scripts )

add_custom_target( geosx_generate_install_schema
         ALL
         COMMAND geosx -s ${CMAKE_BINARY_DIR}/schema.xsd >generate_schema.log 2>&1 || (cat generate_schema.log && exit 1)
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
         DEPENDS geosx
         COMMENT "Generating XML schema to install"
         )

if (NOT GEOSX_INSTALL_SCHEMA)
  set_target_properties(geosx_generate_install_schema PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
endif()

install( FILES ${CMAKE_BINARY_DIR}/schema.xsd
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/schema OPTIONAL)


################################
# Add integratedTests
################################
if( EXISTS "${CMAKE_SOURCE_DIR}/../integratedTests/CMakeLists.txt")
  add_subdirectory( ${CMAKE_SOURCE_DIR}/../integratedTests integratedTests )
else()
  message( "Could not find the integratedTests submodule" )
endif()
