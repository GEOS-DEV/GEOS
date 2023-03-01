

add_custom_target( pygeosx_unit_tests
                   COMMAND ${CMAKE_BINARY_DIR}/lib/PYGEOSX/bin/python3 -m pytest -v
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   DEPENDS pygeosx 
                 )

