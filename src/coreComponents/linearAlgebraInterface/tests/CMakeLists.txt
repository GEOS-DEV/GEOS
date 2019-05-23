#
# Specify list of tests
#

set( LAI_tests
     testLAOperations.cpp 
     testArrayLAOperations.cpp )

set( LAI_tests_no_opts
     testDofManager.cpp )

set( LAI_tests_no_opts_routines
     createConnLocPattern.cpp )

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
#  Add tests without external parameters
#

foreach(test ${LAI_tests_no_opts})
  get_filename_component( test_name ${test} NAME_WE )
  message("processing test_name ${test_name}")
        
  blt_add_executable( NAME ${test_name}  
                      SOURCES ${LAI_tests_no_opts_routines} ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON geosx_core gtest linearAlgebraInterface )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name}
                NUM_MPI_TASKS ${nranks} )

endforeach()
