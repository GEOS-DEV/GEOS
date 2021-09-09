################################
# documentation targets
################################

blt_add_doxygen_target( geosx_doxygen )

# Required directory for doxygen to generate output to
add_custom_target( build_doxygen_output_dir
                   ALL
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/html
                 )

add_dependencies( geosx_doxygen build_doxygen_output_dir )

# Add a CTest test to check for Doxygen warnings
set( mkdir_cmd "${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/html" )
add_test( NAME testDoxygenCheck
          COMMAND bash -c "${mkdir_cmd} && ${DOXYGEN_EXECUTABLE} Doxyfile 2> >(tee doxygen.err) >/dev/null && exit $(cat doxygen.err | wc -l)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          )
