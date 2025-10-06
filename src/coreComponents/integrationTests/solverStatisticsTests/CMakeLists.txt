# Specify list of tests
set( solverStatisticsTests_sources
      testSolverStats.cpp )

set( dependencyList mainInterface  gtest ${parallelDeps} )

# Add gtest C++ based tests
foreach(test ${solverStatisticsTests_sources})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()
