#
# Specify list of tests
#

set( linearAlgebraInterface_tests
     testEpetraLAOperations.cpp )

#
# Add gtest C++ based tests
#

foreach(test ${linearAlgebraInterface_tests})
	get_filename_component( test_name ${test} NAME_WE )
	message("processing test_name ${test_name}")
	blt_add_executable( NAME ${test_name}
	                    SOURCES ${test}
	                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
	                    DEPENDS_ON ${geosx_core_libs} gtest linearAlgebraInterface )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} )

endforeach()
