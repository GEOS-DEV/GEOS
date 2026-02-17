#[[
Package: fluidFlowSolvers unitTests
Builds unit test for fluid flow solver statistics
#]]

set( fluidFlow_gtest_tests
     testVelocityReconstruction.cpp
     testFlowStatistics.cpp )

set( dependencyList mainInterface testingUtilities gtest ${parallelDeps} )

foreach( test_source ${fluidFlow_gtest_tests} )
  get_filename_component( test_name ${test_source} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test_source}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
                      FOLDER UnitTests )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
