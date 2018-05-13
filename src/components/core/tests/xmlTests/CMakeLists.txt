#
# Specify list of tests
#

set(gtest_geosx_tests
   testXML.cpp
   )

#
# Add gtest C++ based tests
#
set(componentlist PhysicsSolverPackage1 geosx_core )
set(componentlistlink "")
foreach(component ${componentlist})
  set( componentlistlink ${componentlistlink} ${GEOSX_LINK_PREPEND_FLAG} ${component} ${GEOSX_LINK_POSTPEND_FLAG} ) 
endforeach()



foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${componentlistlink} gtest ${thirdPartyLibs}  "-lz"
                        )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/basic_input.xml
                  )
endforeach()
