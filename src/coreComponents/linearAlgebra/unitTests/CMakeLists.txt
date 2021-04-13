#
# Specify list of tests
#
set( LAI_tests
     testLAOperations.cpp 
     testArrayLAOperations.cpp
    )

set( nranks 2 )

#
# Add gtest C++ based tests
#
set( dependencyList gtest linearAlgebra )

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

foreach(test ${LAI_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}  
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

if ( ENABLE_MPI )
  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} -x ${nranks}
                NUM_MPI_TASKS ${nranks} )
else()
  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} )
endif()
endforeach()
