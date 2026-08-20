# Unit tests for fileIO/Outputs

set( gtest_outputs_tests
     testMemoryStats.cpp )

if( ENABLE_VTK )
  list( APPEND gtest_outputs_tests testVTKPolyDataWriter.cpp )
endif()

# Link against mainInterface to resolve basicSetup and ProblemManager symbols
set( dependencyList mainInterface gtest ${parallelDeps} )

foreach(test ${gtest_outputs_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
