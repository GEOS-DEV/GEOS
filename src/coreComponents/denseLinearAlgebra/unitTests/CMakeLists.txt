set( serial_tests
     BlasLapack )

set( dependencyList gtest denseLinearAlgebra )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

foreach(test_name ${serial_tests})
  set( exec_name test${test_name} )
  blt_add_executable( NAME ${exec_name}
                      SOURCES ${exec_name}.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  geos_add_test( NAME testDenseLinearAlgebra_${test_name}
                 COMMAND ${exec_name} )
endforeach()
