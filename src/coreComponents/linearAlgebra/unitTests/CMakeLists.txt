set( serial_tests
     testLinearSolverParametersEnums.cpp
     testComponentMask.cpp
     testHypredrive.cpp )

set( parallel_tests
     testMatrices.cpp
     testVectors.cpp
     testExternalSolvers.cpp
     testKrylovSolvers.cpp
     testReverseCutHillMcKeeOrdering.cpp )

set( nranks 2 )

# A single visible HIP device cannot safely service multiple HYPRE MPI ranks:
# HYPRE's ROCm device-side radix sort reports hipErrorIllegalState when ranks
# share the device during matrix assembly. Keep these unit tests runnable on
# the one-GPU HIP configuration while retaining the two-rank CPU coverage.
if( ENABLE_HIP )
  set( nranks 1 )
endif()

set( dependencyList gtest linearAlgebra ${parallelDeps} )
# MPI is intentionally consumed through linearAlgebra's public interface.
# Adding the legacy BLT `mpi` target here independently can cause CMake to
# split libmpicxx and libmpi across the transitive link line.
if( ENABLE_MPI )
  list( REMOVE_ITEM dependencyList mpi )
endif()

foreach(test ${serial_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()

foreach(test ${parallel_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  if ( ENABLE_MPI )
    if( test_name STREQUAL "testExternalSolvers" AND CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
      # SuperLU_DIST's OpenMP direct solve is not stable with AOCC/Clang's
      # runtime in this two-rank unit test. Keep the MPI coverage serial per
      # rank; GEOS OpenMP remains enabled for the rest of the test suite.
      geos_add_test( NAME ${test_name}
                     COMMAND ${test_name} -x ${nranks}
                     NUM_MPI_TASKS ${nranks}
                     NUM_OMP_THREADS 1 )
    else()
      geos_add_test( NAME ${test_name}
                     COMMAND ${test_name} -x ${nranks}
                     NUM_MPI_TASKS ${nranks} )
    endif()
  else()
    geos_add_test( NAME _${test_name}
                   COMMAND ${test_name} )
  endif()
endforeach()
