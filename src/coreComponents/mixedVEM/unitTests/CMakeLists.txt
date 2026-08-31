set( gtest_geosx_tests
     testMixedVEMElementOperators.cpp
     testHybridMixedVEM.cpp )

set( dependencyList mixedVEM gtest ${parallelDeps} )

foreach( test ${gtest_geosx_tests} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
