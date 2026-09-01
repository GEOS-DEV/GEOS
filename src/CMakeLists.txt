
cmake_minimum_required( VERSION 3.24 )

# CUDA architectures are supplied by the host configuration or CI.
cmake_policy(SET CMP0104 NEW)
# BLT 0.6.2 still uses FindCUDA while configuring its CUDA targets.
if(POLICY CMP0146)
  cmake_policy(SET CMP0146 OLD)
endif()
cmake_policy(SET CMP0074 NEW) # dont ignore <PackageName>_ROOT env vars when searching for packages via find_package()
cmake_policy(SET CMP0066 NEW) # use CMAKE_<LANG>_FLAGS_<BUILD_TYPE> for try_compile() instead of only CMAKE_<LANG>_FLAGS
cmake_policy(SET CMP0056 NEW) # use CMAKE_EXE_LINKER_FLAGS in try_compile() in addition to CMAKE_<LANG>_FLAGS

# A legacy geos-tpl installation can provide a generated host config.  The
# filename identifies its compiler and accelerator toolchain. Load it before
# project() so compiler, MPI, and accelerator settings are present when CMake
# enables the corresponding languages.
if( DEFINED GEOS_TPL_DIR AND IS_DIRECTORY "${GEOS_TPL_DIR}" )
  file( GLOB _geos_tpl_named_host_configs CONFIGURE_DEPENDS
        "${GEOS_TPL_DIR}/host-config-*.cmake" )
  list( LENGTH _geos_tpl_named_host_configs _geos_tpl_host_config_count )
  if( _geos_tpl_host_config_count GREATER 1 )
    string( JOIN "\n  " _geos_tpl_host_config_choices ${_geos_tpl_named_host_configs} )
    message( FATAL_ERROR
             "Multiple generated GEOS TPL host configs were found in ${GEOS_TPL_DIR}:\n  ${_geos_tpl_host_config_choices}\nSelect one explicitly with -C <host-config>." )
  elseif( _geos_tpl_host_config_count EQUAL 1 )
    list( GET _geos_tpl_named_host_configs 0 _geos_tpl_host_config )
  elseif( EXISTS "${GEOS_TPL_DIR}/host-config.cmake" )
    # Keep accepting host-config.cmake from older geos-tpl installations.
    set( _geos_tpl_host_config "${GEOS_TPL_DIR}/host-config.cmake" )
  endif()
  if( DEFINED _geos_tpl_host_config AND EXISTS "${_geos_tpl_host_config}" )
    message( STATUS "Loading generated TPL host config: ${_geos_tpl_host_config}" )
    include( "${_geos_tpl_host_config}" )
  endif()
  unset( _geos_tpl_host_config )
  unset( _geos_tpl_named_host_configs )
  unset( _geos_tpl_host_config_count )
  unset( _geos_tpl_host_config_choices )
endif()

################################
# GEOSX
################################
project( geosx LANGUAGES C CXX )
include(GNUInstallDirs)
set( BLT_CXX_STD "c++20" CACHE STRING "Version of C++ standard" )
if(CMAKE_BUILD_TYPE EQUAL "Debug")
  set( ENABLE_WARNINGS_AS_ERRORS "OFF" CACHE PATH "")
else()
  set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "")
endif()
set( CMAKE_ENABLE_EXPORTS ON )

# Configure the arguments used by generated build-system test targets before
# BLT enables CTest. Users who need full control can still define
# CMAKE_CTEST_ARGUMENTS directly from their host-config or CMake command line.
set( GEOS_CTEST_PARALLEL_LEVEL "AUTO" CACHE STRING
     "CTest parallel level for generated test targets. Use AUTO for CTest's default level, OFF for serial execution, or a nonnegative integer." )
set_property( CACHE GEOS_CTEST_PARALLEL_LEVEL PROPERTY STRINGS AUTO OFF 0 1 2 4 8 16 )
if( NOT DEFINED CMAKE_CTEST_ARGUMENTS )
  string( TOUPPER "${GEOS_CTEST_PARALLEL_LEVEL}" _geos_ctest_parallel_level )
  if( _geos_ctest_parallel_level STREQUAL "AUTO" )
    set( CMAKE_CTEST_ARGUMENTS "--parallel" )
  elseif( _geos_ctest_parallel_level STREQUAL "OFF" )
    set( CMAKE_CTEST_ARGUMENTS "" )
  elseif( GEOS_CTEST_PARALLEL_LEVEL MATCHES "^[0-9]+$" )
    set( CMAKE_CTEST_ARGUMENTS "--parallel" "${GEOS_CTEST_PARALLEL_LEVEL}" )
  else()
    message( FATAL_ERROR "GEOS_CTEST_PARALLEL_LEVEL must be AUTO, OFF, or a nonnegative integer." )
  endif()
  unset( _geos_ctest_parallel_level )
endif()

# populated by `geos_add_test` (private list of all geos tests, keep in sync with the macro definition)
set_property( GLOBAL PROPERTY geos_tests_exe_list )

################################
# Include blt
################################
option( ENABLE_BENCHMARKS "Enables benchmarks" ON )
include( cmake/blt/SetupBLT.cmake )
set( BLT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/cmake/blt/ )

# ROCm 7.2 reports an illegal-state error when device-side printf is used on
# the gfx10/gfx11 targets.  BLT's HIP smoke tests intentionally exercise
# device-side printf, so leave those tests available as an opt-in while
# keeping the rest of the HIP smoke coverage usable on these GPUs.
option( GEOS_ENABLE_HIP_DEVICE_PRINTF_TESTS
        "Run BLT HIP smoke tests that use device-side printf" OFF )
if( ENABLE_HIP AND ENABLE_TESTS AND NOT GEOS_ENABLE_HIP_DEVICE_PRINTF_TESTS
    AND CMAKE_HIP_ARCHITECTURES MATCHES "(^|;)gfx(10|11)" )
  set_tests_properties( blt_hip_smoke blt_hip_gtest_smoke
                        DIRECTORY ${BLT_BUILD_DIR}/tests/smoke
                        PROPERTIES DISABLED TRUE )
endif()

# BLT's bundled GoogleTest currently uses an implicit char8_t-to-char32_t
# conversion that Clang 22 diagnoses as -Wcharacter-conversion. Keep GEOS'
# warnings-as-errors policy intact while limiting the compatibility flag to
# the third-party test target.
if( CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 22 )
  foreach( _geos_gtest_target gtest gtest_main gmock gmock_main )
    if( TARGET ${_geos_gtest_target} )
      target_compile_options( ${_geos_gtest_target} PRIVATE -Wno-error=character-conversion )
    endif()
  endforeach()
endif()

# use, i.e. don't skip the full RPATH for the build tree
set( CMAKE_SKIP_BUILD_RPATH FALSE )

# when building, don't use the install RPATH already
# (but later on when installing)
set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
set( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
endif()

# Removes -pthread from BLT CUDA libraries lists
if( ENABLE_CUDA )
  unset( _cuda_libraries )
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
unset( externalComponentsLinkList )
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
message( STATUS "Configuring GEOSX version ${GEOS_VERSION_FULL}" )


################################
# Create header of configuration options
################################
include( cmake/GeosxConfig.cmake )

################################
# Add main
################################
unset( extraComponentsLinkList )

if( ENABLE_OPENMP )
  list( APPEND extraComponentsLinkList openmp )
endif()

if( ENABLE_CUDA )
  list( APPEND extraComponentsLinkList cuda )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND extraComponentsLinkList CUDA::nvToolsExt )
endif()

if( ENABLE_HIP )
  list( APPEND extraComponentsLinkList blt::hip )
endif()

blt_add_executable( NAME geosx
                    SOURCES main/main.cpp
                    DEPENDS_ON geosx_core 
                               ${extraComponentsLinkList}
                               ${externalComponentsLinkList} )

if(APPLE)
  target_link_options(geosx PRIVATE "-Wl,-export_dynamic")
endif()


# Seems to be required on some CMake versions (e.g. 3.16) to get enforce device linking
if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA"  )
  set_target_properties( geosx PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS TRUE )
endif()


# To change the runtime path during installation
set_target_properties( geosx PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
set_target_properties( geosx PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

install(TARGETS geosx RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
# Duplicate the installed `geosx` executable and rename it `geos` to let the end users adapt their scripts.
install( FILES ${CMAKE_INSTALL_PREFIX}/bin/geosx
         DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
         RENAME geos
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

if( ENABLE_XML_UPDATES AND ENABLE_MPI AND UNIX AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )

  set(SCHEMA_DIR ${CMAKE_SOURCE_DIR}/coreComponents/schema)
  set(SCRIPT_DIR ${CMAKE_SOURCE_DIR}/../scripts)
  set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/..")

  add_custom_target( geosx_generate_schema
                     ALL
                     COMMAND geosx -s ${SCHEMA_DIR}/schema.xsd >generate_schema.log 2>&1 || (cat generate_schema.log && exit 1)
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx
                     COMMENT "Generating XML schema"
                     )

  add_custom_target( geosx_build_datastructure_tables
                     ALL
                     COMMAND python3 ${SCRIPT_DIR}/SchemaToRSTDocumentation.py -s ${SCHEMA_DIR}/schema.xsd -o ${CMAKE_SOURCE_DIR}/docs/sphinx >update_rst_tables.log 2>&1 || (cat update_rst_tables.log && exit 1)
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_generate_schema
                     COMMENT "Generating schema-dependent RST files"
                     )

  add_custom_target( geosx_validate_all_xml_files
                     ALL
                     COMMAND bash ${SCRIPT_DIR}/validateXMLFiles.bash ${SCHEMA_DIR}/schema.xsd ${SOURCE_DIR}/examples ${SOURCE_DIR}/inputFiles
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_generate_schema
                     COMMENT "Validating all XML files in the repository against the schema"
                     )
  if( SPHINX_FOUND )
    blt_add_sphinx_target( geosx_docs )
    add_dependencies( geosx_docs geosx_build_datastructure_tables )
  endif()
endif()

if( ENABLE_PYGEOSX )
   add_subdirectory( pygeosx )
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

if( NOT GEOS_INSTALL_SCHEMA )
  set_target_properties(geosx_generate_install_schema PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
endif()

install( FILES ${CMAKE_BINARY_DIR}/schema.xsd
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/schema OPTIONAL)


################################
# Add python environment setup
################################
message(STATUS "Using geosPythonPackages branch feature/remove-dndx-detj")
set(GEOS_PYTHON_PACKAGES_BRANCH "feature/remove-dndx-detj" CACHE STRING "" FORCE)


if ( Python3_EXECUTABLE )
    message(STATUS "Found python version ${Python3_VERSION}")
    if (${Python3_VERSION} VERSION_LESS "3.6.0")
      message(STATUS "Note: try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config to the appropriate version.")
      message( FATAL_ERROR "Building the GEOSX python tools requires Python >= 3.6." )
    endif()

    # Select the version of python to target
    if( ENABLE_PYGEOSX )
        set( PYTHON_POST_EXECUTABLE ${CMAKE_BINARY_DIR}/lib/PYGEOSX/bin/python CACHE PATH "" FORCE )
    else()
        set( PYTHON_POST_EXECUTABLE ${Python3_EXECUTABLE} CACHE PATH "" FORCE )
    endif()

    # Check for the virtualenv package
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import virtualenv"
        RESULT_VARIABLE VIRTUALENV_AVAILABLE
    )

    if (NOT ${VIRTUALENV_AVAILABLE} EQUAL 0)
        message(WARNING "The \"virtualenv\" package was not found in the target python environment (${Python3_EXECUTABLE}).  This package may be required to build PYGEOSX or the python development environment.")
    endif()

    # Check for the requested branch of geosPythonPackages (used for testing)
    if( NOT GEOS_PYTHON_PACKAGES_BRANCH )
      set(GEOS_PYTHON_PACKAGES_BRANCH "main" CACHE STRING "" FORCE)
    endif()

    set( GEOS_PYTHON_TOOLS_BINS
         "${CMAKE_BINARY_DIR}/bin/preprocess_xml"
         "${CMAKE_BINARY_DIR}/bin/format_xml" )

    add_custom_command( OUTPUT ${GEOS_PYTHON_TOOLS_BINS}
                        COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/setupPythonEnvironment.bash -p ${PYTHON_POST_EXECUTABLE} -b ${CMAKE_BINARY_DIR}/bin --python-pkg-branch ${GEOS_PYTHON_PACKAGES_BRANCH}
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

    add_custom_target( geosx_python_tools
                       DEPENDS ${GEOS_PYTHON_TOOLS_BINS} )

    add_custom_target( geosx_python_tools_clean
                       COMMAND rm ${GEOS_PYTHON_TOOLS_BINS} )

    add_custom_target( geosx_python_tools_test
                       COMMAND ${CMAKE_BINARY_DIR}/python/geosx/bin/test_geosx_xml_tools
                       COMMAND rm -r ${CMAKE_BINARY_DIR}/python/geosx_xml_tools_tests*
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
                       DEPENDS geosx_python_tools
                     )

    add_custom_target( geosx_format_all_xml_files
                       COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/formatXMLFiles.bash -g ${CMAKE_BINARY_DIR}/bin/format_xml ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../examples
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS geosx_python_tools
                     )

else()
    message(WARNING "Building the GEOSX python tools requires Python >= 3.6.")
    message(STATUS "If you need these, try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config.")
endif()


################################
# Add integratedTests
################################
if (NOT DEFINED ENABLE_ATS)
    set( ENABLE_ATS true CACHE BOOL "")
endif()


if ( ENABLE_ATS )
    if (NOT DEFINED ATS_WORKING_DIR)
        message( STATUS "ATS_WORKING_DIR is not defined; defaulting to ${CMAKE_BINARY_DIR}/integratedTests/workingDir" )
        set( ATS_WORKING_DIR "${CMAKE_BINARY_DIR}/integratedTests/workingDir" CACHE PATH "")
    endif()

    if (NOT DEFINED ATS_BASELINE_DIR)
        message( STATUS "ATS_BASELINE_DIR is not defined; defaulting to ${CMAKE_SOURCE_DIR}/../integratedTests" )
        set( ATS_BASELINE_DIR "${CMAKE_SOURCE_DIR}/../integratedTests" CACHE PATH "")
    endif()

    if (NOT Python3_EXECUTABLE)
        message( FATAL_ERROR "An appropriate version of python was not found (required for integrated testing system).  Try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config." )
    endif() 

     # Setup testing
    set( ATS_SCRIPT
        "${CMAKE_BINARY_DIR}/integratedTests/geos_ats.sh"
       )

    add_custom_command( OUTPUT ${ATS_SCRIPT}
                        COMMAND ${CMAKE_BINARY_DIR}/bin/setup_ats_environment ${CMAKE_SOURCE_DIR}/.. ${CMAKE_BINARY_DIR} ${ATS_BASELINE_DIR} ${ATS_WORKING_DIR} ${ATS_ARGUMENTS}
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       )

    add_custom_target( ats_environment
                       DEPENDS geosx_python_tools
                       DEPENDS ${ATS_SCRIPT} )

    set( _geos_ats_launch ${ATS_SCRIPT} )
    if( GEOS_ENABLE_SANITIZERS )
      set( _geos_ats_launch ${CMAKE_COMMAND} -E env ${GEOS_SANITIZER_CMAKE_ENV} ${ATS_SCRIPT} )
    endif()

    add_custom_target( ats_run
                       COMMAND ${_geos_ats_launch}
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )

    add_custom_target( ats_clean
                       COMMAND ${_geos_ats_launch} -a veryclean
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )

    add_custom_target( ats_rebaseline
                       COMMAND ${_geos_ats_launch} -a rebaseline
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )

    add_custom_target( ats_rebaseline_failed
                       COMMAND ${_geos_ats_launch} -a rebaselinefailed
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )
    unset( _geos_ats_launch )
endif()


# Python formatting
if ( ENABLE_YAPF )
    set( integrated_tests_python_sources )
    file( GLOB_RECURSE integrated_tests_python_sources "${CMAKE_SOURCE_DIR}/../inputFiles/*.py" )
    set( integrated_tests_ats_sources )
    file( GLOB_RECURSE integrated_tests_ats_sources "${CMAKE_SOURCE_DIR}/../inputFiles/*.ats" )

    blt_add_code_checks( PREFIX  integrated_tests_yapf_style
                         SOURCES ${integrated_tests_python_sources} ${integrated_tests_ats_sources} ${CMAKE_SOURCE_DIR}/coreComponents/dummy.cpp
                         YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )
endif()


# the following adds a `build_test` CMake target such that running `$ make build_test test`
# builds the unit tests before running them
get_property( tmp GLOBAL PROPERTY geos_tests_exe_list )
add_custom_target( build_test COMMAND ${CMAKE_CTEST_COMMAND} ${CMAKE_CTEST_ARGUMENTS} DEPENDS ${tmp} )

geos_apply_sanitizer_test_environment()

message( DEBUG "geos_tests: ${tmp} " )
