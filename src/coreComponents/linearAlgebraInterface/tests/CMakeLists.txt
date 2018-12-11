#
# Specify list of tests
#

set( linearAlgebraInterface_tests
     testDofManager.cpp       
     testLAOperations.cpp )

set( nranks 2 )

#
# Add gtest C++ based tests
#

foreach(test ${linearAlgebraInterface_tests})
	get_filename_component( test_name ${test} NAME_WE )
	message("processing test_name ${test_name}")
	blt_add_executable( NAME ${test_name}  
	                    SOURCES ${test}
	                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
	                    DEPENDS_ON geosx_core gtest linearAlgebraInterface )

  #TODO: be smarter about distinguishing tests that need XML input from those that don't:

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/testDofManager.xml -x ${nranks}
                NUM_MPI_TASKS ${nranks} )

endforeach()
