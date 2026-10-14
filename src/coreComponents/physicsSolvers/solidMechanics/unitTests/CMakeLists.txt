#[[
Package: solidMechanicsSolvers unitTests
Builds unit tests for the solid mechanics solvers
#]]

set( solidMechanics_gtest_tests
     testModalAnalysis.cpp )

set( dependencyList mainInterface gtest ${parallelDeps} )

foreach( test_source ${solidMechanics_gtest_tests} )
  get_filename_component( test_name ${test_source} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test_source}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
                      FOLDER UnitTests )
  # Meshes read by the tests
  target_compile_definitions( ${test_name} PRIVATE
                              "GEOS_MODAL_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\"" )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

  if( ENABLE_MPI )
    set( nranks 2 )
    geos_add_test( NAME ${test_name}_mpi
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  endif()
endforeach()
