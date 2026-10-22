# SPDX-License-Identifier: LGPL-2.1-only
#
# Unit tests for linearAlgebra utilities (LAIHelperFunctions)

set( parallel_tests
     testLAIHelperFunctions.cpp )

set( nranks 2 )

# HYPRE's HIP device assembly is not safe when multiple MPI ranks share one
# visible GPU. Keep the two-rank coverage for CPU builds, but run this test
# with one rank in the single-GPU HIP configuration.
if( ENABLE_HIP )
  set( nranks 1 )
endif()

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
