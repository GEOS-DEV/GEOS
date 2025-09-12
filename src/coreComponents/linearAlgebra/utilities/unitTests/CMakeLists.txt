# SPDX-License-Identifier: LGPL-2.1-only
#
# Unit tests for linearAlgebra utilities (LAIHelperFunctions)

set( parallel_tests
     testLAIHelperFunctions.cpp )

set( nranks 2 )

set( dependencyList mainInterface gtest ${parallelDeps} )

foreach(test ${parallel_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  if ( ENABLE_MPI )
    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  else()
    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )
  endif()
endforeach()
