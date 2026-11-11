#[[
Package: mixedMimetic unitTests
Builds the unit tests of the cell classification (ConsistencyAdaptation)
#]]

set( mixedMimetic_gtest_tests
     testConsistencyAdaptation.cpp )

set( dependencyList mainInterface gtest ${parallelDeps} )

foreach( test_source ${mixedMimetic_gtest_tests} )
  get_filename_component( test_name ${test_source} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test_source}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
                      FOLDER UnitTests )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

  if( ENABLE_MPI AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )
    set( nranks 2 )
    geos_add_test( NAME ${test_name}_mpi
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  endif()
endforeach()
