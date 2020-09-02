

# If an explicit python post-processing target is not set, then try
# using the python being used with cmake.  Note: the target python must
# have the virtualenv package available.
if( NOT PYTHON_POST_EXECUTABLE )
    set( PYTHON_POST_EXECUTABLE ${PYTHON_EXECUTABLE} CACHE PATH "" FORCE )
endif()


set(GEOSX_XML_TOOLS_BINS
  "${CMAKE_BINARY_DIR}/bin/preprocess_xml"
  "${CMAKE_BINARY_DIR}/bin/format_xml"
)


add_custom_command( OUTPUT ${GEOSX_XML_TOOLS_BINS}
                    COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/setupVirtualPythonEnvironment.bash -p ${PYTHON_POST_EXECUTABLE} -o ${CMAKE_BINARY_DIR}/python
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/python/geosx/bin/preprocess_xml ${CMAKE_BINARY_DIR}/bin/preprocess_xml
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/python/geosx/bin/format_xml ${CMAKE_BINARY_DIR}/bin/format_xml
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pygeosx.so ${CMAKE_BINARY_DIR}/python/geosx/lib/python3.8/site-packages/pygeosx.so
                    )

add_custom_target( geosx_xml_tools 
                   DEPENDS ${GEOSX_XML_TOOLS_BINS}
                   )

add_custom_target( geosx_xml_tools_test
                   COMMAND ${CMAKE_BINARY_DIR}/python/geosx/bin/test_geosx_xml_tools
                   COMMAND rm -r ${CMAKE_BINARY_DIR}/python/geosx_xml_tools_tests*
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
                   DEPENDS ${GEOSX_XML_TOOLS_BINS}
                   )

add_custom_target( geosx_format_all_xml_files
                   COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/formatXMLFiles.bash -g ${CMAKE_BINARY_DIR}/bin/format_xml ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../examples
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                   DEPENDS ${GEOSX_XML_TOOLS_BINS}
                   )

