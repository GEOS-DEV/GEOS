# Specify list of tests
set( solverStatisticsTests_sources
      testSolverStats.cpp
      testCoupledSolverStats.cpp )

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

# Check that only one rank writes the statistics CSV files
if( ENABLE_MPI )
    set( nranks 2 )
    geos_add_test( NAME testCoupledSolverStats_mpi
                   COMMAND testCoupledSolverStats
                   NUM_MPI_TASKS ${nranks} )
endif()
