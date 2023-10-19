

if ( Python3_EXECUTABLE )
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
         message(WARNING "To build the PYGEOSX and/or the GEOSX python tools, \"virtualenv\" should be installed in the target python.")
         message(STATUS "To install virtualenv, try running ${PYTHON_POST_EXECUTABLE} -m pip install virtualenv")
     endif()

     # Build targets
     set( GEOSX_PYTHON_TOOLS_BINS
          "${CMAKE_BINARY_DIR}/bin/preprocess_xml"
          "${CMAKE_BINARY_DIR}/bin/format_xml"
        )

     add_custom_command( OUTPUT ${GEOSX_PYTHON_TOOLS_BINS}
                         COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/setupPythonEnvironment.bash -p ${PYTHON_POST_EXECUTABLE} -b ${CMAKE_BINARY_DIR}/bin
                         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       )

     add_custom_target( geosx_python_tools
                        DEPENDS ${GEOSX_PYTHON_TOOLS_BINS} )

     add_custom_target( geosx_python_tools_test
                        COMMAND ${CMAKE_BINARY_DIR}/python/geosx/bin/test_geosx_xml_tools
                        COMMAND rm -r ${CMAKE_BINARY_DIR}/python/geosx_xml_tools_tests*
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
                        DEPENDS geosx_python_tools
                      )

     add_custom_target( geosx_format_all_xml_files
                        COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/formatXMLFiles.bash -g ${CMAKE_BINARY_DIR}/bin/format_xml ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../examples
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                        DEPENDS geosx_xml_tools
                      )

else()
     message(WARNING "Building the GEOSX python tools requires Python >= 3.7.")
     message(STATUS "If you need these, try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config.")
endif()


# Python formatting
if ( ENABLE_YAPF )
     set( python_module_sources )
     file( GLOB_RECURSE python_module_sources "*.py" )

     # Note: blt throws an error if sources doesn't include a c-file, so include dummy.cpp
     blt_add_code_checks( PREFIX  python_modules_yapf_style
                          SOURCES ${python_module_sources} ${CMAKE_SOURCE_DIR}/coreComponents/dummy.cpp
                          YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )

     set( python_script_sources )
     file( GLOB_RECURSE python_script_sources "${CMAKE_SOURCE_DIR}/../scripts/*.py" )

     blt_add_code_checks( PREFIX  python_scripts_yapf_style
                          SOURCES ${python_script_sources} ${CMAKE_SOURCE_DIR}/coreComponents/dummy.cpp
                          YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )
endif()


