#
# Specify list of tests
#

if(ENABLE_PVTPackage)
set( gtest_geosx_tests
     testCompMultiphaseFlow.cpp
     testCompMultiphaseFlowCapPressure.cpp )
endif(ENABLE_PVTPackage)

set( gtest_geosx_tests_with_xml
     testLaplaceFEM.cpp )

set( nranks 2 )

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )
  message("processing test_name ${test_name}")
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON geosx_core gtest )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} ${CMAKE_CURRENT_LIST_DIR} )
endforeach()

#
# Add tests requiring input xml
#
foreach(test ${gtest_geosx_tests_with_xml})
  get_filename_component( test_name ${test} NAME_WE )
  message("processing test_name ${test_name}")
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON geosx_core gtest )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/${test_name}.xml -x ${nranks}
                NUM_MPI_TASKS ${nranks} )
endforeach()
