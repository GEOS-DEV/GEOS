#
# Specify list of tests
#

set( LAI_tests
     testLAOperations.cpp )

set( LAI_tests_with_xml
     testDofManager.cpp )

set( nranks 2 )

#
# Add gtest C++ based tests
#

foreach(test ${LAI_tests})
  get_filename_component( test_name ${test} NAME_WE )
  message("processing test_name ${test_name}")

  blt_add_executable( NAME ${test_name}  
	              SOURCES ${test}
	              OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
	              DEPENDS_ON gtest linearAlgebraInterface )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} -x ${nranks}
                NUM_MPI_TASKS ${nranks} )

endforeach()

# 
#  Add tests requiring input xml
#

foreach(test ${LAI_tests_with_xml})
  get_filename_component( test_name ${test} NAME_WE )
  message("processing test_name ${test_name}")
	
  blt_add_executable( NAME ${test_name}  
                      SOURCES ${test}
	              OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
	              DEPENDS_ON geosx_core gtest linearAlgebraInterface )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/${test_name}.xml -x ${nranks}
                NUM_MPI_TASKS ${nranks} )

endforeach()


