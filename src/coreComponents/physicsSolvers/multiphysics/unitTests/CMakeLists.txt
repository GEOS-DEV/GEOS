set( multiphysics_gtest_tests
       testCoupledSolverStatisticsCSV.cpp )

  set( dependencyList mainInterface testingUtilities gtest ${parallelDeps} )

  blt_add_executable( NAME testCoupledSolverStatisticsCSV
                      SOURCES testCoupledSolverStatisticsCSV.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
                      FOLDER UnitTests )

  geos_add_test( NAME testCoupledSolverStatisticsCSV
                 COMMAND testCoupledSolverStatisticsCSV )

  if( ENABLE_MPI )
    geos_add_test( NAME testCoupledSolverStatisticsCSV_2ranks
                   COMMAND testCoupledSolverStatisticsCSV
                   NUM_MPI_TASKS 2 )
  endif()
