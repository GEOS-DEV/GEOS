#
# Specify list of tests
#

set(gtest_geosx_tests
   testUnstructured.cpp
   )

#
# Add gtest C++ based tests
#
set(componentlist geosx_core )
set(componentlistlink "")
foreach(component ${componentlist})
  set( componentlistlink ${componentlistlink} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} ) 
endforeach()



foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON geosx_core gtest optionparser  "-lz"
                        )

    blt_add_test( NAME ${test_name}_hex 
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/cube_hex.xml
                  )

    blt_add_test( NAME ${test_name}_tet 
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/cube_tet.xml
                  )

    blt_add_test( NAME ${test_name}_pyr 
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/cube_pyr.xml
                  )

    blt_add_test( NAME ${test_name}_wedge
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/cube_wedge.xml
                  )
endforeach()
