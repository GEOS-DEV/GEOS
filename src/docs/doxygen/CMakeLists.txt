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
                   
